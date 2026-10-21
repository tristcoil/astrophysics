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
cd astrophysics
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
| `reboundx`     | General-relativistic and tidal forces for REBOUND (N1, M4)   |
| `gaiamock`     | Gaia epoch astrometry and DR3's non-single-star cascade, for selection functions; GitHub only, needs GSL (M2) |
| `sympy`        | Computer algebra: metrics, Christoffel symbols, Killing equations (N5, N6) |
| `camb`         | CMB power spectra, recombination and neutrino effects (D3, O4, P1, P3) |
| `uproot`, `awkward` | Reading CERN ROOT files in pure Python (P7)             |
| `particle`     | PDG particle masses, spins and quantum numbers (P2, P6)      |
| `aiapy`, `demregpy` | AIA temperature responses and differential emission measure maps (F10) |
| `pyspedas`, `cdasws` | Solar wind and spacecraft data from NASA's SPDF/CDAWeb (F8, F9) |

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
| A5 | Planets being born (ALMA discs)           | Exoplanets     | ★★☆       | ~2.5 GB     | Rings and gaps carved by young planets |
| A6 | Watch planets orbit another star          | Exoplanets     | ★★☆       | ~115 MB     | HR 8799's four planets in motion |
| B1 | Hertzsprung–Russell diagram from Gaia     | Stars          | ★☆☆       | ~50–200 MB  | All of stellar evolution in one plot |
| B2 | 3D map of the solar neighborhood          | Stars          | ★★☆       | ~70 MB      | Fly through real nearby stars |
| B3 | Hunt star clusters by their motion        | Stars          | ★★☆       | ~50 MB      | Clusters pop out of the noise |
| B4 | Cepheids and the Leavitt law              | Stars          | ★☆☆       | ~5–200 MB   | The first step of the distance ladder |
| B5 | Asteroseismology: weigh stars by their sound | Stars          | ★★☆       | ~1.1 GB     | Hear red giants ring |
| B6 | How old are the oldest stars?             | Stars          | ★★☆       | ~2.1 GB     | Can stars be older than the universe? |
| B7 | Reading starlight: spectroscopy           | Stars/galaxies | ★★☆       | ~50–500 MB  | Spectral classes, redshifts, the BPT diagram |
| B8 | A time-lapse of stellar evolution in clusters | Stars      | ★★☆       | ~0.3–1 GB   | Watch stars age across 7,200 clusters |
| B9 | Eclipsing binaries: weigh and measure stars | Stars        | ★★☆       | <50 MB      | Masses and radii to 2% |
| B10 | White dwarfs: the cooling embers         | Stars          | ★★☆       | ~0.1–0.5 GB | Dead stars freezing into crystals |
| B11 | Run a star to core collapse (MESA, simulation) | Stars    | ★★★       | a few GB (code) | The onion inside a star before it explodes |
| B12 | Constellations in motion                  | Stars          | ★☆☆       | ~0.1 GB     | The sky ±100,000 years, and stars that brush past the Sun |
| B13 | Stellar atmospheres from spectra         | Stars          | ★★☆       | ~0.2–1 GB   | Temperature, gravity and Payne's hydrogen |
| B14 | Starspots and cycles on other suns       | Stars          | ★★☆       | ~0.1–1 GB   | Is the Sun unusually quiet? |
| B15 | Stellar coronae in X-rays (eROSITA)      | Stars          | ★★☆       | <100 MB     | 138,800 hot coronae and the activity limit |
| B16 | Stellar winds: Sun to massive stars      | Stars          | ★★☆       | <100 MB     | P Cygni profiles and 3,000 km/s gales |
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
| D5 | Is dark energy changing? (DESI BAO)       | Cosmology      | ★★☆/★★★   | ~0.15 GB/~3.2 GB | The 2025 hint, reproduced |
| D6 | Fresh supernovae onto our Hubble diagram  | Cosmology      | ★★☆       | ~0.4 GB     | Test D1 on supernovae it never saw |
| D7 | Is the universe really expanding?         | Cosmology      | ★☆☆       | ~0.1 GB     | Three tests that tired light fails |
| D8 | The universe's speedometer: H(z)          | Cosmology      | ★☆☆       | <1 MB       | When the expansion began to speed up |
| D9 | H0 from lensed quasar time delays         | Cosmology      | ★★☆       | ~190 MB     | A distance ladder with no rungs |
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
| F5 | Anatomy of a sunspot (Hinode)             | The Sun        | ★★☆       | ~0.5–2 GB   | Why a magnetic field makes a dark spot |
| F6 | The magnetic cycle and polar reversals    | The Sun        | ★★☆       | ~0.2–0.5 GB | The 22-year Hale cycle and the next forecast |
| F7 | Flares, superflares and tree-ring storms  | The Sun        | ★★☆       | <100 MB     | How big can a solar flare get? |
| F8 | The solar wind: Parker to Parker Solar Probe | The Sun     | ★★☆       | ~0.1–1 GB   | From the corona to the edge of interstellar space |
| F9 | CMEs and 150 years of space weather       | The Sun        | ★★☆       | <50 MB      | How often a Carrington storm hits |
| F10 | The million-degree corona                | The Sun        | ★★☆       | ~0.2–1 GB   | Temperature maps of the corona |
| F11 | Helioseismology: the Sun rings           | The Sun        | ★★☆       | ~0.1–0.5 GB | Seeing inside the Sun with sound |
| F12 | Does the Sun's brightness change?        | The Sun        | ★☆☆       | <10 MB      | Spots darken, faculae brighten: net +0.1% |
| G1 | Color composites of real galaxies         | Imaging        | ★★☆       | ~510 MB     | Your own pretty pictures, honestly made |
| G2 | Meteorite landings world map              | Fun/easy       | ★☆☆       | ~10 MB      | Quick, satisfying, geographic |
| H1 | Weigh the black hole at the Galactic Centre | Milky Way      | ★★☆       | <1 MB       | The 2020 Nobel, from S2's orbit |
| H2 | Our own rotation curve from 21-cm hydrogen | Milky Way      | ★★☆       | ~0.1–4 GB   | Our galaxy's rotation and spiral arms |
| H3 | Galactic archaeology with Gaia            | Milky Way      | ★★★       | ~0.2–1 GB   | Stellar streams and the galaxy we ate |
| H4 | The Sun's acceleration from quasars       | Milky Way      | ★★☆       | ~0.6 GB     | Feel the Galaxy's pull in 1.6 million quasars |
| H5 | The Milky Way across the spectrum         | Milky Way      | ★☆☆       | <100 MB     | Ten all-sky maps, radio to gamma rays |
| H6 | Fly through the dust around the Sun       | Milky Way      | ★★☆       | ~8.1 GB     | The Local Bubble and the Radcliffe Wave in 3D |
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
| M1 | Weigh stars with binary orbits (Sirius B) | Binary stars   | ★☆☆       | ~40 MB (+570 MB optional) | Kepler's third law weighs a white dwarf |
| M2 | A census of binary stars (Gaia, SB9)      | Binary stars   | ★★☆       | ~0.7 GB     | Is the lonely Sun typical? |
| M3 | Wide binaries: a test of gravity          | Binary stars   | ★★★       | ~3.6 GB     | Newton vs MOND in a million star pairs |
| M4 | Triple stars and the three-body problem   | Binary stars   | ★★☆       | ~0.1–1 GB   | Three suns, and eclipses that switch off |
| M5 | Stars that touch, trade mass and merge    | Binary stars   | ★★☆       | ~0.75 GB    | A merger watched, a merger mispredicted |
| M6 | White dwarf pairs for LISA                | Binary stars   | ★★☆       | ~0.3 GB     | A 7-minute orbit shrinking by gravitational waves |
| N1 | Precessing orbits: Mercury to S2          | Relativity     | ★★☆       | ~125 MB     | Einstein's 43″ per century, reproduced |
| N2 | Bending starlight: 1919 to Cassini        | Relativity     | ★☆☆       | ~70 MB (+1.3 GB optional) | Redo the eclipse that made Einstein famous |
| N3 | Clocks in gravity: GPS to white dwarfs    | Relativity     | ★★☆       | ~3.4 GB (0.45 GB kept) | Watch a satellite clock feel gravity |
| N4 | Geodetic precession and frame dragging    | Relativity     | ★★☆       | ~0.16 GB (0.9 GB all centres) | Spacetime dragged by the spinning Earth, measured |
| N5 | Falling into a black hole (simulation)    | Relativity     | ★★☆       | ~0.36 GB (0.22 GB kept) | The real sky, seen falling into Sgr A* |
| N6 | Killing vectors and conserved quantities  | Relativity     | ★★☆       | 27 MB       | The symmetries behind orbits and shadows |
| N7 | Is Einstein still right? A scorecard      | Relativity     | ★☆☆/★★★   | ~0.2 GB     | A century of tests on one chart |
| O1 | Build a universe: the Friedmann equations | Expanding universe | ★☆☆   | ~0.14 GB    | Big Crunch, eternal expansion, or no Big Bang |
| O2 | The cosmic light cone and horizons        | Expanding universe | ★☆☆   | ~5 MB       | Galaxies receding faster than light |
| O3 | Why is the night sky dark?                | Expanding universe | ★★☆   | <10 MB      | Olbers' paradox, measured from beyond Pluto |
| O4 | Recombination: the universe turns transparent | Expanding universe | ★★☆ | <1 MB      | Why the CMB was released at 3,000 K |
| O5 | The dark ages and the 21-cm signal        | Expanding universe | ★★☆   | <50 MB      | The first stars' fingerprint, and the EDGES puzzle |
| O6 | The cosmic dipole                         | Expanding universe | ★★☆   | ~2.5 GB     | Is the universe the same in every direction? |
| P1 | Inflation on trial                        | Early universe | ★★☆       | ~0.23 GB    | Which inflation models survive the CMB? |
| P2 | The thermal history: quarks to atoms      | Early universe | ★☆☆       | ~3 MB       | A calendar of the first second |
| P3 | The cosmic neutrino background            | Early universe | ★★☆       | ~0.95 GB    | 336 relic neutrinos per cm³, and their mass |
| P4 | Freeze-out: the WIMP miracle              | Early universe | ★★☆       | ~30 MB      | Why there is as much dark matter as there is |
| P5 | Where did the antimatter go?              | Early universe | ★☆☆       | ~3 MB       | One baryon per 1.6 billion photons survived |
| P6 | The quark–gluon plasma (ALICE)            | Early universe | ★★★       | ~9 MB       | The universe at 20 µs, recreated at the LHC |
| P7 | The Higgs boson in CMS open data          | Early universe | ★★☆       | ~15 GB      | Find the Higgs yourself |

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

**Sixth wave** (sections M, N, O and P, added 2026-10-04): binary and triple stars, general relativity, the Friedmann universe and the Big Bang, and the particle physics of the first second.
- *Binary and triple stars:* M1 → M2 → B9 → M4 → M5 → M6, with M3 as the stretch goal; C9, I4 and I11 are the compact-object versions.
- *Testing Einstein:* N7 (the scorecard) → N2 (light) → N3 (clocks) → N1 (orbits) → N4 (gyroscopes), then C9, C1 and I7 for strong gravity.
- *Inside a black hole's geometry:* N6 (symmetries) → N5 (falling in) → I6 (ray tracing) → I5 (the real shadow) → I7 (the area theorem).
- *The Big Bang, forward in time:* P1 (inflation) → P2 (the thermal history) → P7 (the Higgs) → P6 (quarks bind) → P5 (antimatter) → P3 (neutrinos) → P4 (dark matter freeze-out) → K1 (nuclei) → O4 (atoms) → O5 (the first stars) → J8 (the first galaxies).
- *The expanding universe:* O1 (Friedmann) → O2 (horizons) → D7 (is it expanding?) → D8 (H(z)) → O6 (is it the same in every direction?) → O3 (why the sky is dark).
- *Quick wins (★☆☆):* M1, N2, N7, O1, O2, P2, P5.

**Seventh wave** (F5–F12 and B13–B16, added 2026-10-04): the Sun's surface, atmosphere and wind, and the atmospheres of other stars.
- *The Sun from inside out:* F11 (interior) → F5 (a sunspot) → F6 (the magnetic cycle) → F10 (the corona) → F7 (flares) → F9 (CMEs) → F8 (the wind), with F12 for the Sun's brightness.
- *The Sun as a star:* F12 → B14 (spots and cycles) → B15 (coronae) → F7 (superflares) → B16 (winds) → B13 (atmospheres).
- *Quick win (★☆☆):* F12.

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

**Data.** The EPIC archive and API (`epic.gsfc.nasa.gov`; no key needed): natural-colour and enhanced images with metadata (spacecraft, Sun and Moon positions) for 3,646 days from 13 June 2015 to October 2026 (55,179 images), typically 13–22 images a day. The calibrated ten-band images (level 1B, ~385 MB each) are archived at NASA's Atmospheric Science Data Center behind a free NASA Earthdata Login; there is no public mirror. The A4 build used only the public images (see `44-earth-as-exoplanet/README.md`).

**Method + visualizations.**
- The keeper: a year of Earth in one movie, with an "Earth as a dot" panel beside it showing the dot's colour and brightness.
- Light curves in each band; autocorrelation and periodogram → the rotation period. How few images per day would still recover it?
- Colour against time: is the Pacific-facing hemisphere darker and bluer than the Africa-facing one? How much does cloud cover change the answer?
- Seasons: the North Pole tips into view in June, the South Pole in December.
- The Moon transit of 16 July 2015 as an animation, and the Moon's brightness against Earth's (the Moon reflects ~12% of sunlight, Earth ~30%).

**Check against.** Jiang et al. 2018 (AJ 156, 26): the 24-hour period recovered by autocorrelation.

**Watch out for.** The Sun–Earth–spacecraft angle swings between 1.8° and 12.4° (2015–2026, measured from the API's vectors; the often-quoted ~4°–15° is too high) on an 87.5-day cycle as DSCOVR loops around L1. The dot brightens by ~0.9% per degree toward small angles, ~4× a Lambert sphere; correct for it. The natural-colour images are processed for display with an unpublished, frame-to-frame-varying tone curve, so use level-1B data for photometry. The archive has gaps: a safe mode from 27 June 2019 to 12 February 2020, and 15 July to 9 October 2025.

**Builds on.** A1 and A2 (exoplanets as points of light), E2 (orbits), F2 (another planet-sized disk, imaged every day).

## A5. Planets being born: ALMA's protoplanetary discs ★★☆

**Science question.** Young stars are surrounded by discs of gas and dust from which planets form. ALMA's sharpest images show these discs carved into rings, gaps and spirals. Are we seeing young planets clearing their paths, and how massive would they have to be?

**Background.** A planet embedded in a disc pushes gas away from its orbit and opens a gap, and dust grains collect in rings where the gas pressure peaks. A gap's width and depth grow with the planet's mass relative to the star, and depend on the disc's viscosity and temperature. In 2014 ALMA's first long-baseline image showed concentric rings around HL Tau, a star about a million years old (ALMA Partnership 2015). DSHARP (Disk Substructures at High Angular Resolution Project; Andrews et al. 2018) imaged 20 nearby discs at 35 mas (~5 au) resolution in 1.25 mm continuum and CO. Rings and gaps are everywhere; a few discs have spirals. From the gap widths, Zhang et al. (2018) inferred planets from about Neptune's to Jupiter's mass.

**Data.** The DSHARP data release (`almascience.eso.org/almadata/lp/DSHARP/`): fiducial continuum images of all 20 discs as FITS (35–62 MB each, 0.8 GB in all), DSHARP II's own radial profiles, CO line cubes (0.5–1.9 GB each; HD 163296's is 1.7 GB) and the calibrated measurement sets. HL Tau's ALMA long-baseline images (ALMA science verification data, 16 MB). Gaia distances to each star (Gaia DR3 has no parallax for HL Tau, and Elias 27's is too poor: use their clouds' distances).

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

**Background.** Direct imaging suppresses the star's light (with a coronagraph and image processing) to reveal planets tens of thousands of times fainter right beside it. HR 8799, 41 pc away and ~40 Myr old, has four planets of several Jupiter masses at ~15–70 au. The outermost takes about 500 years per orbit and the innermost about 55. Over a decade their motion is plain to see. Measuring each planet's position relative to the star (astrometry, to ~1 mas with the Gemini Planet Imager) constrains Keplerian orbits. Wang et al. (2018) found that the system is stable only if the four orbits are nearly coplanar and close to a 1:2:4:8 chain of period resonances. β Pictoris b, on a ~20-year orbit seen nearly edge-on, is a second example; its data ship with the `orbitize!` package.

**Data.** Published astrometry of the four HR 8799 planets compiled by Wang et al. (2018, AJ 156, 192; GPI points with ~1 mas errors), plus earlier Keck measurements (e.g. Konopacky et al. 2016). β Pic b astrometry from `orbitize!`'s example data. Gaia DR3 for the star's distance, and the Hipparcos–Gaia Catalog of Accelerations (Brandt 2021) for a dynamical mass constraint.

**Method + visualizations.**
- Fit each planet's orbit with `orbitize!` (the OFTI method for short arcs, MCMC otherwise) → posterior orbits, periods and inclinations. (As built: OFTI fails on the 1-mas GPI points of d and e, so the main fit is a joint nested-sampling fit of all four planets, with `orbitize!` as the cross-check.)
- The keeper: an animation of the four planets moving along their fitted orbits from 1998 to 2050, with the measured positions as dots.
- Test coplanarity and the 1:2:4:8 chain; run posterior draws through N-body integrations (REBOUND, as in E1) — which survive 40 Myr?
- β Pic b's edge-on orbit, and its closest passages in front of the star.

**Check against.** Wang et al. 2018: orbit fits and stability results.

**Watch out for.** Short arcs leave orbital elements strongly correlated: quote posteriors, not single best fits. Different instruments have different plate scales and north angles; fit offsets between them (Wang et al.'s GPI points need De Rosa et al. 2020's +0.36° north correction). Planet masses come from evolutionary models and the age unless the dynamics constrain them. The arXiv source archives the tables are parsed from total ~115 MB (1 MB of LaTeX is kept).

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

**Data.** The OGLE-IV Collection of Variable Stars, classical Cepheids of the Magellanic System (Soszyński et al. 2015, extended to 2023: 9,667 stars, LMC 4,713 + SMC 4,954). It provides periods, mean I and V magnitudes, pulsation modes, and full I/V light curves, via `https://www.astrouw.edu.pl/ogle/ogle4/OCVS/` (`lmc/cep`, `smc/cep`; the same tree as the FTP server). Milky Way Cepheids with Gaia DR3 parallaxes (`gaiadr3.vari_cepheid` joined to `gaia_source`) set the zero point; OGLE's vetted list of Galactic Cepheids (`OCVS/allGalCep.listID`, Pietrukowicz et al. 2021, updated 2026) gives a clean sample with Gaia DR3 IDs, and C-MetaLL IX (VizieR J/A+A/708/A216) and Trentin et al. 2024 (J/A+A/681/A65) give the spectroscopic metallicities the zero point needs. Leavitt's original 25 stars (Harvard College Observatory Circular 173, 1912) are public domain; the scan is at `https://articles.adsabs.harvard.edu/pdf/1912HarCi.173....1L`, and OGLE's identification file lists their Harvard (HV) numbers.

**Method + visualizations.**
- A few folded light curves (the classic sawtooth rise and slow decline).
- Leavitt's 25 stars next to OGLE's thousands: the same law, nearly 400× the sample.
- Period–luminosity in V, I and W for the LMC. Fundamental-mode and first-overtone pulsators form two parallel sequences. The scatter drops from ~0.2 mag in V to ~0.08 mag in W.
- SMC vs LMC offset → relative distance. The SMC's larger scatter shows its depth along the line of sight.
- Residuals from the law mapped across the LMC show which side of its tilted disk is nearer.
- Zero point from Gaia parallaxes of Milky Way Cepheids → absolute LMC distance.

**Check against.** LMC distance 49.59 ± 0.54 kpc (distance modulus 18.477; Pietrzyński et al. 2019, from eclipsing binaries, a method independent of ours). Wesenheit slope ≈ −3.3 (Soszyński et al. 2015).

**Watch out for.** Separate the pulsation modes before fitting. Gaia parallaxes have a zero-point offset (Lindegren et al. 2021 recipe; bright Cepheids may need more, as Riess et al. 2021 found). Milky Way Cepheids are more metal-rich than LMC ones, so state the metallicity term we assume. Blending in crowded fields makes Cepheids look too bright. (Found when building it: the metallicity term moves the LMC distance by ~0.17 mag, so it has to be fitted from spectroscopic [Fe/H]; the parallax-space fit with a free parallax offset is biased by regression dilution, here by 0.036 mag, so test it on simulated calibrators; Gaia's Wesenheit coefficient 1.90 over-corrects heavily reddened Milky Way Cepheids; and Gaia's catalogue-mean G of bright Cepheids is 0.05–0.08 mag fainter than their light-curve intensity mean, so use `vari_cepheid`'s `int_average_*`.)

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

**Watch out for.** The scaling relations need percent-level corrections (APOKASC-3 applies them; say whether we do). Long cadence cannot see νmax above ~250 µHz. Frequency resolution needs months of continuous data. Granulation noise dominates at low frequency. The solar reference values depend on the pipeline: B5's code measures νmax☉ = 3143 µHz on SOHO/VIRGO, 1.7% above the SYD pipeline's 3090 µHz (6% in every mass), so calibrate on the Sun with the same code or adopt one catalogue's whole scale. A plain autocorrelation biases the Δν of red-clump stars by ~0.5%.

**Builds on.** A1 (same archive and tools), B1 (the red clump, seen from the inside).

**Extensions.** Red clump vs red giant branch from mixed-mode period spacings (Bedding et al. 2011): stars of identical surface appearance, burning helium or not. Solar-like oscillations of main-sequence stars in Kepler short cadence (e.g. 16 Cyg A/B).

## B6. How old are the oldest stars? ★★☆

**Science question.** Globular clusters hold the Milky Way's oldest stars. How old are they? And are they younger than the universe that D1 and D3 measured?

**Background.** A cluster's stars formed together. As a cluster ages, the **main-sequence turnoff** (where stars exhaust core hydrogen) moves to fainter, redder stars. Fitting theoretical isochrones to the turnoff reads off the age. In the 1990s globular cluster ages (~15 Gyr) seemed older than a universe with a high H0 allowed: the "age crisis". It was resolved by Hipparcos distances and dark energy. Our own projects restage it. D1's supernovae give 12.5 Gyr with H0 = 73.5; D3's CMB fit gives 13.825 ± 0.033 Gyr.

**Data.** The HST UV Legacy Survey of Galactic Globular Clusters (HUGS; Piotto et al. 2015, Nardiello et al. 2018): astro-photometric catalogues of 56 globular clusters in five bands (F275W, F336W, F438W, F606W, F814W) with membership probabilities, a MAST High Level Science Product. MIST isochrones with HST bolometric corrections: use MIST v2.5 ("MIST II", Dotter et al. 2026), which adds α-enhanced models ([α/Fe] up to +0.6), rather than the scaled-solar v1.2 (Choi et al. 2016). Its ACS/WFC package is one 1.1 GB tarball, `data/tarballs_v2.5/isos/HST_ACS_WFC.txz` at mist.science (the site's `HST_ACSWF.txz` link is broken), and it lacks [Fe/H] = −2.25 and −2.75. Distances from Baumgardt & Vasiliev (2021, Gaia + HST; the per-method table is in the arXiv source, not on VizieR); metallicity and reddening from the Harris catalogue (2010 edition).

**Method + visualizations.**
- Clean colour–magnitude diagrams (quality flags, membership) of M92, NGC 6397, M15 and 47 Tuc.
- Fix distance and reddening from independent sources, then fit the turnoff and subgiant branch on an age grid. Monte Carlo over distance, reddening, metallicity and [α/Fe] → an age distribution for each cluster.
- Age vs metallicity: the metal-rich 47 Tuc is younger than the metal-poor clusters.
- The ultraviolet colours split each cluster into multiple populations (the "chromosome map"). The survey was designed for this.
- The ages on a timeline against D1's and D3's ages of the universe.

**Check against.** M92: 13.80 ± 0.75 Gyr (Ying et al. 2023), at their own distance modulus of 14.72 ± 0.04. Baumgardt & Vasiliev's 14.65 makes M92 ~0.8 Gyr older.

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

*B13–B16 look at the outer layers of other stars: their atmospheres, spots, hot coronae and winds. The Sun (F5–F12) is the one star whose atmosphere we can resolve in detail; these projects ask how typical it is.*

## B13. Stellar atmospheres: temperature, gravity and composition from spectra ★★☆

**Science question.** A star's spectrum forms in a thin atmosphere a few hundred to a few thousand kilometres deep. Model atmospheres compute how light escapes through it, and fitting them to a spectrum gives the star's temperature, surface gravity and composition. In 1925 Cecilia Payne used the first such physics to show that stars are made mostly of hydrogen, against the consensus of the day. Can we fit model atmospheres to real spectra across the HR diagram, tell giants from dwarfs by their line shapes, and test the models' predicted limb darkening against exoplanet transits?

**Background.** In a stellar atmosphere temperature falls outward, so absorption lines form where the gas is cooler than the layers below. The Saha and Boltzmann equations set which ions and energy levels are populated; B7 used them to explain why Balmer lines peak near 9,500 K. Payne's thesis showed that the huge differences between stellar spectra come from temperature, not composition, and that hydrogen and helium vastly outnumber everything else; Russell, who had first dismissed the result, confirmed it in 1929. Surface gravity sets the gas pressure: in a dwarf (log g ≈ 4.5) collisions broaden strong lines into wide wings, while a giant (log g ≈ 2) has narrow lines and stronger lines of ionized metals. Grids of model atmospheres and synthetic spectra, such as PHOENIX (Husser et al. 2013), cover 2,300–12,000 K. The same models predict limb darkening, how a star's disk dims toward its edge, which shapes every exoplanet transit (A1) and was measured for the Sun in F2.

**Data.** The X-shooter Spectral Library DR3 (Verro et al. 2022): 830 flux-calibrated spectra of 683 stars from 350 to 2,480 nm at medium resolution (`xsl.u-strasbg.fr`). The PHOENIX synthetic spectra and specific-intensity grids (`phoenix.astro.physik.uni-goettingen.de`; download only the needed temperatures). Gaia DR3 astrophysical parameters for comparison (I11's queries). Kepler transit light curves from A1, and published limb-darkening coefficients (Claret 2017 for TESS, Claret & Bloemen 2011 for Kepler).

**Method + visualizations.**
- Fit PHOENIX spectra to XSL stars by χ² over a grid in temperature, log g and [Fe/H], with a smooth continuum correction → parameters for each star, against the literature and Gaia DR3.
- Giants against dwarfs: Mg b and Na D wings and the Ca II triplet at the same temperature, with model spectra at log g = 2 and 4.5.
- Payne's argument rebuilt: the observed strength of a hydrogen line and of a calcium line across temperature, against Saha–Boltzmann predictions for different hydrogen-to-metal ratios → only a hydrogen-dominated composition fits.
- Inside a model atmosphere: temperature and pressure against optical depth for a 3,500 K, 5,800 K and 10,000 K star, and where the core and wings of a line form.
- Limb darkening: model intensity profiles against Kepler transit shapes from A1 (fit the transit with the coefficients free) and against F2's solar measurement.

**Check against.** XSL's published parameters (Arentsen et al. 2019); PHOENIX's solar spectrum against F4's Kurucz atlas; Claret's limb-darkening tables.

**Watch out for.** Temperature, gravity and metallicity are correlated in fits, especially at medium resolution; fix one where it is known (log g from Gaia distances and radii). One-dimensional models in local thermodynamic equilibrium miss convection and non-LTE effects that matter at the 0.1 dex level and for strong lines. M dwarfs are full of molecular bands that models reproduce imperfectly. Interstellar reddening changes the continuum shape.

**Builds on.** B7 (spectral classes), F4 (the solar spectrum), F2 (solar limb darkening), A1 (transits), B1 (the HR diagram).

## B14. Starspots and magnetic cycles on other suns ★★☆

**Science question.** Do other stars have spots and 11-year cycles like the Sun? Kepler measured tens of thousands of stars dimming as spots rotate into view, and a fifty-year survey at Mount Wilson recorded the magnetic glow of their chromospheres year by year. Can we measure stellar rotation and activity ourselves, find cycles in other stars, use spin as a clock (gyrochronology), and ask whether the Sun is unusually quiet?

**Background.** Spots rotating across a star make its brightness wobble with the rotation period. Young stars spin fast and are very active; magnetized winds brake them over time, so period grows roughly as √age (Skumanich 1972). Calibrated in clusters, the relation turns a rotation period into an age (Barnes 2007), though spin-down stalls for a while in some stars (Curtis et al. 2020). McQuillan, Mazeh & Aigrain (2014) measured 34,030 rotation periods from Kepler. Magnetic heating of the chromosphere fills in the cores of the Ca II H and K lines; the Mount Wilson S index measures that emission. Wilson's survey, begun in 1966, found Sun-like cycles in many stars, flat activity in others (possibly in Maunder-like minima), and erratic variation in young stars (Baliunas et al. 1995). Compared with Kepler stars of the same temperature and rotation period, the Sun's brightness varies about five times less (Reinhold et al. 2020), so either the Sun is unusually quiet or it can become much more active.

**Data.** McQuillan et al. 2014's rotation periods (VizieR `J/ApJS/211/24`); Kepler light curves via `lightkurve` (A1). Gaia DR3 rotational modulation results (`vari_rotation_modulation`). The Mount Wilson HK Project archive of S-index measurements for ~2,300 stars, 1966–1995, described on the NSO's historical-data pages (the old ftp://solis.nso.edu path is empty; check NSO for the current location), with the Boro Saikia et al. 2018 catalogue of 4,454 cool stars' activity (VizieR) as a fallback. F1's solar cycles and the Sun's S index (Egeland et al. 2017).

**Method + visualizations.**
- Rotation periods ourselves: Lomb–Scargle and autocorrelation for a few hundred Kepler stars → compare with McQuillan's; spot evolution seen as changing light-curve shapes.
- Period against colour for all 34,030 stars: the slow and fast sequences, the gap between them, and lines of constant gyrochronological age; the Kepler open cluster NGC 6811 as a calibrator.
- Mount Wilson time series: S index over decades for stars with clean cycles (cycle lengths 7–14 years), flat stars and erratic young ones, with the Sun's cycle drawn on the same scale.
- Cycle length against rotation period: "active" and "inactive" branches (Böhm-Vitense 2007), and where the Sun falls.
- The Sun among its peers: Kepler variability amplitude against rotation period, with the Sun's amplitude from F12's irradiance data, as in Reinhold et al. 2020.

**Check against.** McQuillan et al. 2014 periods; Baliunas et al. 1995's cycle classes; Reinhold et al. 2020.

**Watch out for.** Kepler's processing removes slow trends, which suppresses periods longer than ~20–30 days and biases amplitudes. Spots at different latitudes rotate at different rates (F1's differential rotation), so measured periods spread. S indices from different instruments need calibration onto one scale. Gyrochronology is less reliable for old stars and for stars whose spin-down stalls.

**Builds on.** F1 and F6 (the solar cycle), F12 (the Sun's variability), A1 (Kepler light curves), B8 (cluster ages), F7 (stellar superflares).

## B15. Stellar coronae: X-rays and the rotation–activity relation ★★☆

**Science question.** The Sun's million-degree corona (F10) glows faintly in X-rays. Other cool stars have coronae too, some a thousand times brighter relative to their light. The eROSITA X-ray telescope found 138,800 coronal stars in its first all-sky survey. Can we map X-ray brightness across the HR diagram, show that it follows rotation, and find the point where a star's magnetism saturates?

**Background.** Stars with outer convection zones, from late F types to M dwarfs, make magnetic fields by dynamo action and heat coronae to 1–10 MK. Hot stars without convective envelopes have no such coronae (their X-rays come from shocks in their winds, B16). The ratio of X-ray to bolometric luminosity, L_X/L_bol, rises steeply with faster rotation, best expressed through the Rossby number Ro = P_rot/τ_c, the rotation period divided by the convective turnover time. Below Ro ≈ 0.1 activity saturates at L_X/L_bol ≈ 10⁻³ (Wright et al. 2011, 824 stars). The Sun, with L_X/L_bol of 10⁻⁷ to 10⁻⁶ over its cycle, sits at the slow, faint end. Fully convective M dwarfs, which lack the Sun's tachocline (F11), follow the same relation (Wright et al. 2016). Because rotation slows with age (B14), X-ray brightness also falls with age, and young planets around active stars receive far more high-energy radiation.

**Data.** The eRASS1 coronal star identifications from the HamStar method (Freund et al. 2024; eROSITA-DE DR1 catalogue pages at `erosita.mpe.mpg.de/dr1/`), matched to Gaia DR3. The Wright et al. 2011 sample with rotation periods (VizieR `J/ApJ/743/48`). B14's Kepler periods for an overlap sample. GOES X-ray fluxes (F2) for the Sun as a star.

**Method + visualizations.**
- The X-ray HR diagram: Gaia colour–magnitude diagram of eRASS1 coronal stars, coloured by L_X/L_bol, showing where convective envelopes switch on (near early-F types) and the young, active populations.
- The rotation–activity relation: L_X/L_bol against Rossby number for Wright's sample plus eRASS1 stars with Kepler or literature periods → the saturation level and the slope of the unsaturated branch.
- The Sun on the same plot, from GOES and published solar soft X-ray luminosities over a cycle.
- Activity against age: X-ray luminosity for stars in open clusters of known age (B8), and the high-energy history of a Sun-like star and of an M dwarf like Proxima.
- Stellar flares in eROSITA light curves or the eRASS1 variability catalogue.

**Check against.** Wright et al. 2011 (saturation at L_X/L_bol ≈ 10⁻³ for Ro ≲ 0.13; unsaturated slope ≈ −2.7); Freund et al. 2024's counts and distributions.

**Watch out for.** eRASS1 covers half the sky and is flux-limited, so active stars dominate at large distances; compare within volume-limited samples. Convective turnover times are model-dependent and set the Rossby scale. Unresolved binaries both spin up stars (tides) and blend X-ray sources. Flares inflate single-epoch X-ray luminosities.

**Builds on.** F10 (the solar corona), B14 (rotation and age), F11 (the tachocline), B8 (cluster ages), A2 (planets around active stars).

## B16. Stellar winds: from the solar breeze to the gales of massive stars ★★☆

**Science question.** The Sun loses about 2 × 10⁻¹⁴ of its mass each year to its wind. Hot massive stars blow winds at up to 3,000 km/s that strip away a large fraction of their mass, and cool giants shed theirs in slow, dusty outflows. Can we measure wind speeds from P Cygni line profiles in archival ultraviolet spectra, see how mass loss scales with luminosity, and detect the faint winds of Sun-like stars through the bubbles they blow in interstellar gas?

**Background.** In a hot star, radiation pushes on the gas through millions of absorption lines, driving a wind (Castor, Abbott & Klein 1975). Its strong ultraviolet lines take a **P Cygni** profile: emission from the whole wind plus blueshifted absorption from gas in front of the star, whose blue edge marks the terminal velocity v∞, typically ~2.6 times the surface escape speed for O and early-B stars (Lamers et al. 1995). Mass-loss rates reach 10⁻⁷–10⁻⁵ M☉/yr, enough to remove much of a massive star's envelope and shape what explodes (L2, B11). Cool giants lose mass through slow winds driven by pulsation and dust; Betelgeuse's Great Dimming of 2019–2020 was a dust cloud formed from a burst of mass loss (Montargès et al. 2021). Sun-like winds are too thin to see directly, but where they meet interstellar gas they build up neutral hydrogen walls, an astrosphere, that absorb a little Lyman-α light. From these, Wood et al. (2005, 2021) found that more active stars have stronger winds, up to an activity limit.

**Data.** IUE high-resolution ultraviolet spectra of O, B and Wolf–Rayet stars, and of cool giants and supergiants, from MAST (`astroquery.mast`, small files). Published terminal velocities (Prinja, Barlow & Howarth 1990) and mass-loss compilations; the astrospheric mass-loss rates of Wood et al. (2021). F8's measured solar wind for the Sun.

**Method + visualizations.**
- P Cygni gallery: C IV, Si IV and N V profiles of O and B stars from IUE, from supergiants to dwarfs, with the blue absorption edge marked.
- Measure v∞ from the edge for ~50 hot stars; compare with Prinja et al., then with escape speeds from published masses and radii → the v∞/v_esc ratio and its jump near 21,000 K (the "bistability jump").
- A simple line-profile model (Sobolev approximation, β velocity law) fitted to one star's C IV doublet.
- Mass loss against luminosity across the HR diagram: O stars, Wolf–Rayet stars, red supergiants, giants and the Sun, with the Reimers (1975) law for cool giants.
- Winds of Sun-like stars: Wood et al.'s mass-loss rates against X-ray surface flux (B15), with the Sun from F8, and what this implies for the young Sun's wind and for planets of active stars.

**Check against.** Prinja et al. 1990 terminal velocities; v∞/v_esc ≈ 2.6 above the bistability jump (Lamers et al. 1995); the Sun's ~2 × 10⁻¹⁴ M☉/yr.

**Watch out for.** Saturated absorption troughs make the blue edge easy to see but turbulence in the wind smears it, which biases v∞ high if ignored. Mass-loss rates from different diagnostics (ultraviolet lines, Hα, radio) disagree by factors of several because winds are clumpy. IUE fluxes need the archive's final calibration (NEWSIPS). Astrospheric absorption depends on interstellar conditions around each star, a large systematic.

**Builds on.** F8 (the solar wind), B15 (stellar activity), B11 (mass loss in stellar evolution), L2 (what explodes), H6 (interstellar gas around the Sun).

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
- A compressed CMB prior on (θ*, Ω_b h², (Ω_b + Ω_c) h²) (Lemos & Lewis 2023; DESI DR2 appendix A), either from D3's own fit or the one DESI used. *(Corrected 2026-10-04: the third parameter is ω_b + ω_c. DESI's released mean for it, 0.14208, includes the 0.06 eV neutrino; the chains give 0.14143 for ω_b + ω_c. See D5's README.)*
- Optional: DES-Y5 supernovae (public on GitHub) to test the sample dependence.
- For our own BAO measurement: BOSS DR12 CMASS/LOWZ large-scale-structure catalogues with randoms (SDSS SAS), or DESI DR1 LSS catalogues (`data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/`).

**Method + visualizations.**
1. BAO alone in ΛCDM → Ωm and r_d·h. DESI's Ωm (≈0.30) vs D1's flat Ωm (0.332 ± 0.018): the mild disagreement that drives the result.
2. Add the CMB prior, then fit w0–wa (MCMC with `emcee`) for BAO + CMB + Pantheon+. Draw the w0–wa contours for each probe and for the combination. Compute Δχ² against ΛCDM → the significance.
3. Swap in DES-Y5 and watch the significance move.
4. Stretch: the correlation function ξ(s) of BOSS CMASS or DESI DR1 LRGs (Landy–Szalay estimator with randoms and weights, reusing D2's pair-count code) → the BAO bump near 100 h⁻¹ Mpc (Eisenstein et al. 2005) → our own D_V/r_d.
- Also: a Hubble diagram of the BAO distances with the ΛCDM and best w0wa curves, and the implied w(z) band.

**Check against.** DESI DR2 ΛCDM Ωm = 0.2975 ± 0.0086; Δχ² = 4.9 (1.7σ) for DESI + Pantheon+ and 8.0 (2.4σ) for DESI + compressed CMB. *(Note 2026-10-04: the 2.8σ for DESI + CMB + Pantheon+ uses the full CMB likelihood with lensing; with the compressed prior, DESI's own released best fits give 2.3σ. DESI's best fits and chains are public at data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/, an exact answer key.)*

**Watch out for.** w0–wa is a convenient parametrization, not a physical model. An apparent "phantom crossing" (w < −1 in the past) may be an artefact of it. The significance depends on supernova calibration (see the debate over DES-Y5's low-redshift sample, e.g. Efstathiou 2025, MNRAS 538, 875). In late 2025 the DES-Dovekie recalibration (Popovic et al. 2025, arXiv:2511.07517) cut DES's 4.2σ to 3.2σ. Quote the result as "Δχ² = X for 2 extra parameters", not as a discovery.

**Builds on.** D1 (likelihood), D2 (clustering code), D3 (CMB prior, sound horizon). Together with B4 and C4, a fourth independent angle on the Hubble tension.

## D6. Fresh supernovae onto our Hubble diagram ★★☆

**Science question.** Pantheon+ was frozen in 2022. Can we fit light curves of Type Ia supernovae that it never saw, and check whether they land on D1's Hubble diagram?

**Background.** Type Ia supernovae are *standardizable*, not standard. Brighter ones fade more slowly (larger "stretch" x1), and redder ones look dimmer (colour c). The SALT light-curve model fits (peak brightness, x1, c) to multi-band photometry. The Tripp relation `μ = m_B − M + α x1 − β c` then gives the distance.

**Data.** The ZTF Type Ia supernova Data Release 2 (Rigault et al. 2025): ~3,600 spectroscopically classified SNe Ia at z < 0.3 from 2018–2020, with public g/r/i light curves, a sample not in Pantheon+. For the live element: recent classifications from the Transient Name Server (TNS; free account for API use), the ZTF Bright Transient Survey explorer, and ZTF photometry from the ALeRCE broker API (public). *(Built 2026-10-04: the DR2 light curves are 130 MB read by byte range out of a 1.37 GB zip, plus 134 MB of SFD dust maps, so ~0.4 GB in all. TNS's public search CSV needs no account; the BTS explorer has classified nothing since August 2026; ALeRCE's `v2/lightcurve/forced-photometry` endpoint serves the forced photometry now inside ZTF alerts, but only about half of recent TNS SNe Ia are in the public stream.)*

**Method + visualizations.**
- SALT3 fits with `sncosmo` (Milky Way dust from the SFD map), with fit quality cuts.
- Standardize with D1's α, β and M → distance moduli → overplot on D1's Hubble diagram and residual plot. *(D1 never fitted α, β itself: they are recovered from Pantheon+'s table, together with its simulated bias corrections. M only applies to fits made with Pantheon+'s own SALT2 model, minus its 0.27 mag `MAG_OFFSET`, which SNANA applies and sncosmo does not. SALT3 + D1's numbers without that puts every supernova 0.19 mag too faint.)*
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
- The DES 5-year supernova release (GitHub `des-science/DES-SN5YR`, folder `0_DATA`): photometry of all 19,706 DES transients and the SALT3 fits of the 1,635 cosmology supernovae at 0.06 < z < 1.12, of which White et al. (2024) kept 1,504 to measure time dilation. Use tag `1.3`: the main branch has since moved to the 2026 "Dovekie" recalibration (photometry unchanged).
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

**Data.** 32 chronometer measurements of H(z) at 0.07 ≤ z ≤ 1.965, compiled in Moresco et al. (2022, *Living Reviews in Relativity*). Moresco's `CCcovariance` repository on GitLab holds precise values for his own 15 and the redshift-dependent systematic errors from which the covariance of all 32 is built (the review applies them to every point; only that reproduces its fit). For comparison: D5's DESI D_H/r_d converted to H(z), and D1's fit.

**Method + visualizations.**
- H(z) against ΛCDM (D1/D3 parameters), a matter-only universe and an empty one.
- ȧ(z) and its dip. A smooth fit or a Gaussian process gives the transition redshift and its uncertainty.
- The deceleration parameter `q(z) = (1 + z) H′/H − 1`, and H0 from the chronometers alone.
- The future: a(t) from the Big Bang to +100 Gyr for ΛCDM (exponential growth), for D5's w0–wa fit, and for w < −1. With w < −1 a "Big Rip" tears everything apart in finite time (Caldwell, Kamionkowski & Weinberg 2003).

**Check against.** D1's transition redshift of 0.59; D5's DESI H(z) points.

**Watch out for.** Stellar-population systematics (initial mass function, metallicity, star-formation history) dominate the error budget, so a diagonal-only fit overstates the precision. The published recipe scales that shared systematic by each measured H value, which pulls fits low (Peelle's puzzle; H0 by ~2 km/s/Mpc here); scale it by a smooth fit instead. Progenitor bias: galaxies selected as passive at different redshifts must belong to the same population. Extrapolating w0–wa into the future is speculation; label it as such.

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

**Background.** The **halo mass function** (Press & Schechter 1974; Tinker et al. 2008) gives the number of dark matter halos per volume above mass M at redshift z. Its high-mass end is exponentially sensitive to Ωm and to σ8, the typical density fluctuation on 8 h⁻¹ Mpc scales. Surveys usually quote `S8 = σ8 (Ωm / 0.3)^0.5`. Several low-redshift measurements came out 2–3σ below the CMB's S8: the **S8 tension**. eROSITA's first all-sky survey (eRASS1) counted 5,259 clusters over 12,791 deg² and found S8 = 0.86 ± 0.01, in line with the CMB (Ghirardini et al. 2024). As of 2026, KiDS-Legacy cosmic shear agrees with the CMB, while DES Y6 sits 2.0σ (cosmic shear, S8 = 0.798) to 2.6σ (3×2pt, S8 = 0.789 ± 0.012) below the combined Planck + ACT + SPT CMB.

**Data.** The eRASS1 cluster catalogues (Bulbul et al. 2024): the primary catalogue of 12,247 optically confirmed clusters and groups (48 MB) and the cosmology subsample, on the eROSITA-DE DR1 site (`erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/`; the plain `.fits` links answer 403, the `.fits.tgz` bundles work). The selection function is *not* released as a file: D11 rebuilt it from the curves plotted by Ghirardini et al. (2024, Fig. 4, a vector PDF) and Clerc et al. (2024), the footprint from Bulbul et al.'s Fig. 2, and the depth from the eRASS1 exposure-map HiPS. Planck's PSZ2 catalogue (1,653 SZ detections) for comparison. Halo mass functions from `colossus`.

**Method + visualizations.**
- Sky map of the eRASS1 clusters (the western Galactic hemisphere) and their redshift distribution.
- Counts in bins of X-ray count rate (the selection observable; not luminosity) and redshift against predictions on an (Ωm, σ8) grid, using the published mass–observable relation and a selection function rebuilt from the published curves.
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

*F5–F12 go below and above the visible surface. The Sun's light comes from the **photosphere**, a layer a few hundred kilometres thick at ~5,800 K, churned by convection into granules. Magnetic fields poking through it make sunspots. Above it the temperature dips to ~4,000 K, then climbs through the **chromosphere** and a thin **transition region** to the million-degree **corona**, which boils off into space as the **solar wind**. Magnetic energy stored in the corona is released in **flares** and **coronal mass ejections**. Suggested solar arc: F5 (a sunspot up close) → F6 (the magnetic cycle) → F11 (inside the Sun) → F10 (the corona) → F7 (flares) → F9 (CMEs) → F8 (the wind) → F12 (the Sun's brightness). B13–B16 compare the Sun with other stars.*

## F5. Anatomy of a sunspot ★★☆

**Science question.** Why are sunspots dark? Up close a sunspot is a dark core (the umbra) ringed by a filamentary penumbra, set in a carpet of boiling granules. In 1908 George Ellery Hale split the light of sunspots into several lines and showed they were magnetic, the first magnetic field found beyond the Earth. Can we measure a sunspot's magnetic field, temperature and flows ourselves from space-telescope data, and see why a strong field makes a cool, dark spot?

**Background.** In a magnetic field each spectral line splits into components with opposite circular polarization (the Zeeman effect), separated in proportion to the field strength. Umbral fields are typically 2,000–3,500 G, thousands of times the Earth's; the strongest measured in a sunspot is 6,250 G (Okamoto & Sakurai 2018). The field suppresses the convection that carries heat up from below, so the umbra is ~4,000 K against the surroundings' 5,800 K. By Stefan–Boltzmann it emits only about a quarter as much light per area, and looks dark only by contrast. Magnetic pressure supports the spot, so its gas is thinner and more transparent and we see deeper into it: the Wilson depression, several hundred kilometres deep, noticed by Alexander Wilson in 1769 from the way spots look foreshortened near the limb. Gas streams outward along the penumbral filaments at a few km/s (the Evershed effect, 1909). Around the spot, granules about 1,000 km across live for 5–10 minutes.

**Data.** Hinode's Solar Optical Telescope: Spectro-Polarimeter Level 2 maps of field strength, inclination, azimuth and Doppler velocity from the MERLIN inversions (JAXA DARTS, `darts.isas.jaxa.jp/solar/hinode/data/sotsp_level2.html`; also the HAO/CSAC archive), and filter images in the G band and Ca II H at 0.1–0.2″ resolution (Hinode Science Data Centre Europe, `sdc.uio.no`). HMI continuum images and magnetograms for the whole disk (F2's routes).

**Method + visualizations.**
- The keeper: a large sunspot in Hinode's G band at full resolution, beside maps of field strength, inclination and line-of-sight velocity on the same grid.
- Profiles along a cut from the umbra to the quiet Sun: field strength and inclination (vertical in the umbra, nearly horizontal at the penumbra's edge), continuum intensity → brightness temperature via the Planck function.
- Field strength against umbral intensity for many spots (darker means stronger), with Hale's 1908 measurement for scale.
- The Evershed flow: Doppler velocities along penumbral filaments, corrected for viewing angle.
- Granulation: cell sizes and lifetimes from a Hinode time series, by segmentation and tracking → the convective velocities that the spot's field suppresses.
- The Wilson depression from the asymmetric shape of a spot as rotation carries it to the limb (HMI continuum over several days).

**Check against.** Umbral fields of 2,000–3,500 G and umbra–photosphere temperature differences of 1,500–2,000 K (Solanki 2003 review); Evershed speeds of a few km/s; granule sizes ~1 Mm and lifetimes 5–10 min.

**Watch out for.** Inversions assume a model atmosphere and give field values averaged over the resolution element, with a filling factor; check the fit quality maps. The azimuth has a 180° ambiguity. Brightness temperatures depend on wavelength and on scattered light in the telescope, which brightens umbrae. Near the limb, foreshortening and the formation height of the lines change what is measured.

**Builds on.** F1 (sunspot records), F2 (HMI magnetograms, limb darkening), F4 (spectral lines), N3 (spectral line shifts).

## F6. The Sun's magnetic cycle: polar reversals and the 22-year Hale cycle ★★☆

**Science question.** The 11-year sunspot cycle is half of a 22-year magnetic cycle: at each solar maximum the Sun's north and south magnetic poles swap. Can we watch the field of decaying sunspots drift to the poles and flip them, measure the tilts of sunspot pairs that drive the flip, and use the polar field to forecast the next cycle?

**Background.** Hale's law: in each hemisphere the leading spots of bipolar regions share one polarity, opposite in the other hemisphere, and the polarities swap every cycle (F2 checked 10 of 10 Cycle 25 regions). Joy's law: the pairs are tilted, the leading spot nearer the equator, more so at higher latitudes. As regions decay, their trailing flux drifts poleward with the meridional flow, cancels the old polar field and builds a new one of opposite sign: the Babcock–Leighton mechanism (Babcock 1961; Leighton 1969). The polar field peaks near solar minimum, and its strength then is the best predictor of the next cycle's amplitude (Schatten et al. 1978; Muñoz-Jaramillo et al. 2013). Wilcox Solar Observatory has measured the polar fields since 1976. Solar Orbiter began viewing the poles from out of the ecliptic in 2025.

**Data.** Wilcox Solar Observatory polar field measurements, 1976 to now (`wso.stanford.edu/Polar.html`). HMI synoptic Carrington-rotation magnetograms, CR 2096 (2010) onward, as FITS files from JSOC's open data directory (`jsoc1.stanford.edu/data/hmi/synoptic/`, no registration). Sunspot-group tilt angles from the Debrecen Photoheliographic Data (`fenyi.solarobs.epss.hu`, DPD). F1's SILSO cycle amplitudes.

**Method + visualizations.**
- The magnetic butterfly diagram: each synoptic map averaged over longitude → field against latitude and time, showing the active-region bands and the poleward surges that reverse the poles.
- WSO's polar fields against time with the sunspot number: the reversals near each maximum, and the 22-year Hale cycle.
- Joy's law from DPD tilts: tilt against latitude, cycle by cycle, and whether stronger cycles have smaller tilts (as some studies find).
- A surface flux transport model (a 1D advection–diffusion equation in latitude), fed with the observed regions → simulated polar field against WSO's. (Labelled simulation.)
- The precursor forecast: polar field at minimum against the next cycle's peak for Cycles 21–25, with Cycle 26's forecast from the current polar field.

**Check against.** Polar reversal times reported for Cycle 24 and 25; the WSO precursor correlation; Hale's and Joy's laws.

**Watch out for.** Polar fields are viewed almost edge-on from the ecliptic, so measurements near the poles are noisy and swing with the Earth's ±7° tilt to the solar equator over the year; WSO's filtered values remove that. A cycle forecast from five or six cycles is statistically weak; show the scatter. HMI synoptic maps fill the far side with older data.

**Builds on.** F1 (butterfly diagram), F2 (Hale's law, magnetograms), F5 (sunspot fields).

## F7. Solar flares, stellar superflares and storms recorded in tree rings ★★☆

**Science question.** How big can a solar flare get, and how often? Fifty years of X-ray monitoring give the statistics for ordinary flares. Thousands of Sun-like stars watched by Kepler show what happens over much longer times. Tree rings record a few storms far bigger than any seen with instruments. Do all three point to the same answer?

**Background.** A flare is the sudden release of magnetic energy stored in twisted coronal fields, through magnetic reconnection. It heats plasma to 10–20 MK, which shines in soft X-rays (F2 caught a 10.6 MK C2 flare). Flares are classed by peak soft X-ray flux measured by GOES: A, B, C, M, X, each ten times the last. Their numbers follow a power law in peak flux, dN/dF ∝ F^(−2) or so (Veronig et al. 2002), so big flares are rare but not separate. The largest instrumentally recorded flare, on 4 November 2003, saturated the detectors; its peak is estimated at about X28 on the old scale, and some estimates go higher. The largest of the current cycle so far was an X9.0 on 3 October 2024. The Carrington event of 1859 was probably larger still. Kepler watched 56,450 Sun-like stars and found 2,889 superflares above 10³⁴ erg, implying such a flare on a Sun-like star about once per century; their energy distribution continues the Sun's (Vasilyev et al. 2024). Tree rings record sudden jumps in radiocarbon made by solar energetic particles: in 774–775 CE, 993–994 CE and at a few earlier dates ("Miyake events", Miyake et al. 2012), the largest 14,300 years ago (Bard et al. 2023).

**Data.** GOES X-ray flare reports 1975–2017 (NOAA NCEI, yearly text files) and the GOES-16+ flare summaries since 2017 (NCEI, F2's routes). The Vasilyev et al. 2024 superflare catalogue (*Science* data supplement) and Kepler light curves for a few of the stars via `lightkurve`. Published annual radiocarbon series across the 774 and 993 CE events (Miyake et al. 2012; Büntgen et al. 2018) and the 14,300-year-old event (Bard et al. 2023).

**Method + visualizations.**
- Flare counts per year by class against F1's sunspot number, 1975–2026.
- The frequency distribution of peak flux → the power-law index, cycle by cycle, with the largest events labelled.
- Flare energy: convert peak X-ray flux to total radiated energy with published scaling relations, so that solar flares and Kepler superflares share an axis.
- The keeper: flare frequency against energy from the Sun (GOES) and from Sun-like stars (Kepler), with the Carrington event and the Miyake events placed by their estimated energies.
- A superflare in a Kepler light curve, fitted with a rise and exponential decay → its energy from the star's luminosity.
- The radiocarbon spikes: Δ¹⁴C against year across 774 CE, with a simple carbon-cycle box model → the size of the particle burst relative to the strongest instrumentally measured events.

**Check against.** Veronig et al. 2002; Vasilyev et al. 2024 (2,889 superflares; once per century above 10³⁴ erg); Miyake et al. 2012.

**Watch out for.** GOES classes changed scale: NOAA scaled GOES-8 to GOES-15 fluxes down by a factor of 0.7 (XRS-B) to match older satellites, and dropped that factor with GOES-16 from late 2019. Science-quality data in true physical units are ~43% higher than the old classes, so put all eras on one scale before counting. Saturated flares need the published corrected peaks. Superflare samples suffer from contamination by close binaries and unresolved neighbours in Kepler's large pixels; Vasilyev et al. used Gaia to vet them. Radiocarbon jumps need a carbon-cycle model to turn into particle numbers, and the energy spectrum of the particles matters.

**Builds on.** F1 (cycles), F2 (GOES, flare temperatures), F3 (a storm in detail), A1 (Kepler light curves), B14 (activity of other stars).

## F8. The solar wind: from Parker's prediction to Parker Solar Probe ★★☆

**Science question.** In 1958 Eugene Parker argued that the Sun's hot corona cannot be held down by gravity and must flow outward faster than sound: a solar wind. Spacecraft confirmed it within four years. Can we solve Parker's equation, test it against six decades of measurements, see the fast wind that pours from the poles, follow Parker Solar Probe into the corona, and find where the wind ends at the edge of interstellar space?

**Background.** For an isothermal corona at temperature T, Parker's equation has one solution that starts slow and passes the sound speed at the critical radius r_c = GM☉m_p/(4kT), about 6 R☉ for 1 MK, reaching hundreds of km/s by the Earth. The Sun's rotation winds the outflowing magnetic field into a spiral, ~45° to the radial direction at 1 AU. At 1 AU the wind averages ~400 km/s, ~5 protons per cm³ and a 5 nT field; it carries away over a million tonnes per second. Ulysses flew over the Sun's poles three times (1994–2008) and found steady ~750 km/s wind from the polar coronal holes at solar minimum, and slow, variable wind at low latitudes (McComas et al. 2000). Parker Solar Probe crossed into the magnetically dominated corona in April 2021 (Kasper et al. 2021) and on 24 December 2024 passed 6.2 million km from the surface at 692,000 km/h. The wind slows abruptly at the termination shock (Voyager 1 in 2004 at 94 AU) and meets the interstellar medium at the heliopause, crossed by Voyager 1 in August 2012 at 122 AU and Voyager 2 in November 2018 at 119 AU.

**Data.** NASA SPDF (via CDAWeb, `cdasws` or `pyspedas`, no registration): hourly OMNI solar wind at the Earth since 1963; Ulysses hourly plasma and magnetic field; Parker Solar Probe SWEAP and FIELDS data; Voyager 1 and 2 magnetic field and Voyager 2 plasma. AIA 193 Å images for coronal holes (F2's routes).

**Method + visualizations.**
- Solve Parker's isothermal equation for 1–2 MK (labelled calculation): speed against distance, with OMNI's 1 AU distribution and Parker Solar Probe's measurements close in.
- Ulysses' famous polar plot: wind speed against latitude around the orbit, at solar minimum and at maximum.
- The Parker spiral: the measured field angle at 1 AU against the angle predicted from the wind speed.
- Coronal holes and fast streams: dark regions in AIA 193 Å images as they face the Earth, against wind speed at Earth 2–4 days later; the 27-day recurrence.
- Six decades of OMNI: speed, density and field against the solar cycle.
- Parker Solar Probe's closest passes: magnetic "switchbacks" (sudden field reversals) and the slowdown of the wind toward the Sun.
- The edge of the heliosphere: Voyager 2's plasma speed dropping at the termination shock, and both Voyagers' magnetic field and particle counts at the heliopause.

**Check against.** McComas et al. 2000 (fast wind ~750 km/s); Kasper et al. 2021 (Alfvén surface crossing at ~19 R☉); the Voyager heliopause crossings (Stone et al. 2013; Richardson et al. 2019).

**Watch out for.** An isothermal wind is a teaching model: the real corona is not isothermal, and the fast wind needs extra heating far from the Sun. OMNI is time-shifted to the bow shock and merged from several spacecraft. Voyager 1's plasma instrument failed in 1980, so its heliopause crossing is seen in particles and plasma waves, not directly in speed. Parker Solar Probe data versions change; record which ones were used.

**Builds on.** F2 (AIA, coronal holes), F3 (OMNI, one storm), F10 (the corona's heat), E10 (comet ion tails point along the wind).

## F9. Coronal mass ejections and 150 years of space weather ★★☆

**Science question.** Several times a day at solar maximum, the Sun throws off a billion-tonne cloud of magnetized plasma. A few hit the Earth and drive geomagnetic storms. How fast do these coronal mass ejections travel, how well can we predict when they arrive, and how often has the Earth been hit by storms like the Carrington event?

**Background.** Coronal mass ejections (CMEs) are seen in coronagraphs, which block the solar disk; SOHO's LASCO has watched since 1996, and the CDAW catalogue lists tens of thousands of CMEs with speeds from ~20 to over 3,000 km/s. "Halo" CMEs appear to surround the Sun because they are heading toward or away from us. In the solar wind a CME slows or speeds up toward the wind's speed by aerodynamic drag (the drag-based model, Vršnak et al. 2013), so fast CMEs arrive sooner than their initial speed alone predicts but later than a constant-speed guess. When a CME's field points southward it couples to the Earth's magnetosphere and drives a storm, measured by the Dst index (since 1957) and the longer aa index (since 1868). The Carrington storm of 1859 arrived 17.6 hours after the flare. In July 2012 a CME at least as fast as Carrington's crossed the Earth's orbit where the Earth had been about a week earlier (Baker et al. 2013).

**Data.** The SOHO/LASCO CDAW CME catalogue (`cdaw.gsfc.nasa.gov/CME_list/`, 1996 onward). The Richardson & Cane list of interplanetary CMEs at the Earth (`izw1.caltech.edu/ACE/ASC/DATA/level3/icmetable2.htm`, updated to 2026). Dst from the Kyoto World Data Center; the aa index from ISGI (`isgi.unistra.fr`). OMNI solar wind (F8).

**Method + visualizations.**
- CME statistics: rates against the solar cycle, distributions of speed, angular width and mass, and halo CMEs.
- Transit times: CDAW speed against travel time to the Earth for CME–ICME pairs, with a constant-speed line, an empirical fit (Gopalswamy et al. 2001) and the drag-based model → which predicts arrival best, and how large the typical error is (about half a day).
- 150 years of storms: the aa index since 1868 and Dst since 1957, with the largest storms ranked (1859 estimates, 1921, 1989, 2003, May 2024).
- How often is "extreme"? The distribution of storm peaks → a recurrence time for Carrington-size storms, with the spread between published estimates.
- The near miss: the July 2012 CME's speed and STEREO-A's in-situ measurements, and its predicted Dst had it hit the Earth.

**Check against.** Gopalswamy et al.'s CDAW statistics; Vršnak et al. 2013; Riley 2012 (~12% chance of a Carrington-size storm per decade, later revised downward by others).

**Watch out for.** Coronagraph speeds are projected on the sky and underestimate the speed of Earth-directed CMEs. Matching CMEs to arrivals is ambiguous when several erupt in a few days (F3). Old geomagnetic indices are not homogeneous (stations, methods), and Carrington's Dst is an estimate from one magnetometer in Bombay. Rare-event extrapolation is very sensitive to the fitted tail; show alternatives.

**Builds on.** F3 (one storm end to end), F7 (flares), F8 (the solar wind).

## F10. The million-degree corona ★★☆

**Science question.** The Sun's surface is 5,800 K, yet the corona above it is over a million degrees. How do we know, and how hot is each part of it? Can we turn the six extreme-ultraviolet channels of SDO into temperature maps of the corona, test the simple physics of coronal loops, and size up the energy problem that has puzzled solar physicists for eighty years?

**Background.** During eclipses the corona shows emission lines that matched no known element ("coronium") until Grotrian (1939) and Edlén (1943) identified them as iron stripped of 9 to 13 electrons, which requires temperatures above a million kelvin. Each AIA extreme-ultraviolet channel is dominated by a few ion species and so responds to a band of temperatures (171 Å near 0.8 MK, 193 Å near 1.5 MK, 335 Å near 2.5 MK, 94 and 131 Å up to 10 MK in flares). Combining six channels gives the **differential emission measure** (DEM): how much plasma lies along each line of sight at each temperature. Active regions reach 3–4 MK, the quiet corona 1–1.5 MK, coronal holes less. In a loop in equilibrium, heating balances conduction and radiation, which links peak temperature, pressure p and length L: T_max ≈ 1,400 (pL)^(1/3) K in cgs units (Rosner, Tucker & Vaiana 1978). The corona needs ~10⁷ erg cm⁻² s⁻¹ of heating over active regions and ~3 × 10⁵ in the quiet Sun (Withbroe & Noyes 1977). Candidates are waves and many small reconnection events ("nanoflares", Parker 1988); Solar Orbiter has seen tiny brightenings dubbed "campfires" (Berghmans et al. 2021).

**Data.** AIA images in the six coronal channels (F2's no-registration routes), with temperature responses from `aiapy`. The `demregpy` regularized DEM inversion (Hannah & Kontar 2012), pure Python. IRIS slit spectra and slit-jaw images of the chromosphere and transition region (Mg II h & k, C II, Si IV) from LMSAL (`iris.lmsal.com/search`). Published eclipse images of the white-light and green-line corona for comparison.

**Method + visualizations.**
- AIA temperature responses, and why each channel "sees" a different layer.
- The keeper: DEM-weighted temperature and emission-measure maps of the full disk, with active regions, quiet Sun and coronal holes.
- Loop scaling: measure lengths and DEM temperatures of active-region loops, densities from emission measure → test the Rosner–Tucker–Vaiana law.
- The atmosphere's temperature from top to bottom: F2's photospheric temperatures, IRIS chromospheric lines and AIA's coronal DEM on one temperature-against-height plot, with a standard model atmosphere (e.g. Avrett & Loeser 2008).
- The energy budget: radiative and conductive losses from the DEM → the heating rate per area, against Withbroe & Noyes.
- Cooling loops after a flare (F2 saw channels peak in temperature order) → cooling times compared with conduction and radiation.

**Check against.** Quiet-Sun ~1–1.5 MK and active-region 2–4 MK DEM peaks; Rosner, Tucker & Vaiana 1978; Withbroe & Noyes 1977.

**Watch out for.** DEM inversion is ill-posed: smooth solutions hide multi-thermal structure, and AIA cannot see plasma below ~0.5 MK or well above ~10 MK. Temperature responses depend on assumed element abundances and atomic data (CHIANTI versions); state them. The corona is optically thin, so every pixel sums plasma along the line of sight, including overlapping loops.

**Builds on.** F2 (AIA, coronal scale heights), F4 (spectral lines), F8 (the wind it drives), F7 (flares).

## F11. Helioseismology: the Sun rings like a bell ★★☆

**Science question.** The Sun's surface rises and falls by a few hundred metres in waves with periods near five minutes: sound waves trapped inside it, ringing in millions of modes. Their frequencies reveal the Sun's interior. Can we find the modes in decades of Sun-as-a-star velocity data, measure the depth of the convection zone, and see how the Sun rotates inside?

**Background.** Sound waves are generated by convection near the surface and travel inward until rising temperature refracts them back; each mode is labelled by its degree ℓ (number of nodal lines on the surface) and radial order n. Low-degree modes, seen in the integrated light of the whole Sun, form a comb of peaks spaced by the large separation Δν ≈ 135 µHz, centred near 3,090 µHz (the same pattern B5 used for red giants). Inverting thousands of mode frequencies gives the sound speed throughout the Sun and the base of the convection zone at 0.713 R☉ (Christensen-Dalsgaard et al. 1991; Basu & Antia 1997). Rotational splitting of the modes shows that the convection zone keeps the surface's differential rotation all the way down, while the radiative interior rotates almost rigidly; the shear layer between them, the tachocline, is a prime suspect for the dynamo (F6). Lowering the solar heavy-element abundance (Asplund et al. 2009) spoiled the agreement between models and helioseismology (the "solar abundance problem"); higher revised abundances (Magg et al. 2022) restore much of it.

**Data.** BiSON (Birmingham Solar Oscillations Network) open data portal (`bison.ph.bham.ac.uk/opendata`): calibrated Sun-as-a-star velocity residuals from 1976 (continuous from the 1990s) and machine-readable published mode frequencies. Model S (Christensen-Dalsgaard et al. 1996), the standard reference solar model, from `users-phys.au.dk/jcd/solar_models/`. Published sound-speed and rotation inversions from MDI and HMI (Schou et al. 1998; Larson & Schou 2018).

**Method + visualizations.**
- The power spectrum of a BiSON year: the comb of p modes around 3 mHz, the large and small separations, and an échelle diagram, as B5 did for Kepler red giants.
- Mode frequencies against the asymptotic relation → Δν and the Sun's mean density, compared with what M☉ and R☉ give.
- Frequency shifts with the solar cycle: the same modes are a few tenths of a µHz higher at solar maximum, tracking F1's sunspot number.
- The interior from published inversions: the relative sound-speed difference between the Sun and Model S against radius, with the convection-zone base; the rotation rate against radius and latitude, showing the tachocline.
- Ray paths: how modes of different ℓ sample different depths, computed from Model S's sound speed.

**Check against.** Δν ≈ 135 µHz, ν_max ≈ 3,090 µHz; base of the convection zone 0.713 ± 0.001 R☉ (Basu & Antia 1997); sound-speed differences below ~0.5% between the Sun and Model S.

**Watch out for.** Daily gaps in ground-based data put sidelobes at ±11.6 µHz around every peak; BiSON's multi-site network reduces them, and the duty cycle matters for fitting. Mode frequencies are fitted with Lorentzian profiles that are slightly asymmetric. The inversions are published results, not ours; full inversions need thousands of modes from resolved Doppler images.

**Builds on.** B5 (asteroseismology of red giants), F1 (the cycle), F6 (the dynamo), K6 and C7 (solar neutrinos test the same solar models).

## F12. Does the Sun's brightness change? ★☆☆

**Science question.** The Sun's total output, about 1,361 W/m² at the Earth, was long called the "solar constant". Space radiometers since 1978 show it is not: it rises and falls with the sunspot cycle by about 0.1%, and a big sunspot group can dim it by a few tenths of a percent for days. Why does a more active Sun shine brighter when sunspots are dark? And how much could solar changes matter for climate?

**Background.** Sunspots block light, but at solar maximum bright magnetic faculae and network, spread over a much larger area, more than make up for them. Models that add facular brightening and subtract sunspot darkening reproduce the measured total solar irradiance (TSI) closely (NRLTSI2, Coddington et al. 2016; SATIRE, Krivova et al.). The absolute level was revised in 2011 to 1,360.8 ± 0.5 W/m² (Kopp & Lean 2011) after scattered light in older radiometers was found. The ultraviolet varies far more than the total: tens of percent in the far ultraviolet, which heats and changes the upper atmosphere. The change in climate forcing over a cycle is about (1 − albedo)/4 × ΔTSI ≈ 0.2 W/m², and the IPCC's best estimate of solar forcing since 1750 is 0.01 W/m², against 2.7 W/m² from human activity (IPCC AR6). Reconstructions of the Maunder Minimum's TSI differ by several times.

**Data.** LASP's LISIRD (`lasp.colorado.edu/lisird`, no registration): TSIS-1 and SORCE TSI measurements, the observational TSI composite, NRLTSI2 with its separate sunspot-darkening and facular-brightening terms, historical reconstructions, and the composite Lyman-α and Mg II indices. F1's sunspot data and F2's HMI images.

**Method + visualizations.**
- TSI since 1978 from each instrument, before and after the 2011 revision, and the composite against F1's sunspot number.
- A spot-and-facula model of our own: sunspot darkening from sunspot areas (F1) and facular brightening from a facular proxy (the Mg II index), fitted to TSIS-1 → how much each contributes.
- A sunspot-group transit: the dip in TSI as a big group crosses the disk (e.g. October 2003 or May 2024), matched with HMI images of the group.
- Spectral variability: the cycle amplitude against wavelength, from visible to Lyman-α.
- Climate in perspective: the solar forcing over a cycle and since 1750, against greenhouse forcing; the range of Maunder Minimum reconstructions.

**Check against.** TSI 1,360.8 ± 0.5 W/m² at solar minimum (Kopp & Lean 2011); a cycle amplitude of ~0.1%; NRLTSI2's components.

**Watch out for.** Instruments disagree on absolute level and drift; composites differ in how they bridge gaps and degradation (the ACRIM–PMOD dispute over a trend between minima). Area-based spot models depend on the sunspot contrast assumed. Climate response to forcing is a separate, model-dependent question; this project quantifies the forcing only.

**Builds on.** F1 (sunspots), F2 (HMI images), F5 (sunspot temperatures), F6 (magnetic cycle).

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

**Data.** Gaia DR3's quasar-like sources of the Gaia-CRF3 (`gaiadr3.agn_cross_id`, ~1.6 million), with their astrometry from `gaia_source`. R₀ from H1. H2 adopted rather than measured V₀ (Reid et al. 2019), so V₀ comes from Sgr A*'s proper motion (Reid & Brunthaler 2020) with H1's R₀. Useful extras: `gaiadr3.frame_rotator_source` (the 428,034 quasars that fix the frame's spin) and Quaia's redshifts (Zenodo 10403370, 171 MB). Download ~0.6 GB with 44 columns per quasar.

**Method + visualizations.**
- Quasar proper motions should be zero: their distribution, centred on zero, with widths set by Gaia's errors.
- Fit vector spherical harmonics (the l = 1 "glide" is the acceleration; the l = 1 rotation is a check on the reference frame), with error weights, outlier rejection and a magnitude- and colour-dependent error model.
- The keeper: an all-sky map of binned proper-motion vectors converging on the Galactic Centre.
- Our acceleration against V₀²/R₀ from H1 and H2, and its direction against the Galactic Centre's.

**Check against.** Klioner et al. 2021: (2.32 ± 0.16) × 10⁻¹⁰ m/s², i.e. 5.05 ± 0.35 µas/yr. Gaia-CRF3 (Klioner et al. 2022): a spin-only fit accepts 1,607,289 sources.

**Watch out for.** The signal is ~90 times smaller than a single quasar's proper-motion error, so only the average over a million sources reveals it. Gaia's own large-scale patterns (frame spin, degree 2–9 swirls and flows) have an rms of ~6 µas/yr, more than the signal: a dipole-only fit points 12° off the centre. Some "quasars" are stars: reject sources with significant parallax or proper motion. Conditioning proper motions on zero parallax needs Gaia's parallax zero point (with zero point 0 the glide jumps to 8 µas/yr), and selecting quasars by parallax biases their proper motions. The drift differs between z < 1.25 and z > 1.25 at 4σ (Quaia redshifts), a systematic of unknown origin. ESA's anonymous TAP queries now stop after ~15 s (sync) or 60 s (async), and async results stop at 50,000 rows unless MAXREC is set.

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

**Background.** Dust dims and reddens starlight. Knowing the distance and the extinction of many stars along each line of sight places the dust in 3D. Edenhofer et al. (2024) built such a map out to 1.25 kpc: Nside 256 (~14′) on the sky, 516 distance bins, parsec-scale in distance. The Sun sits inside the **Local Bubble**, a cavity a few hundred parsecs across blown out by supernovae over the past ~14 Myr. Most young stars nearby formed on its surface (Zucker et al. 2022). Many nearby clouds (Orion, Perseus, Taurus, Cepheus…) line up along the **Radcliffe Wave**, 2.7 kpc long with an amplitude of ~160 pc (Alves et al. 2020), which is also oscillating (Konietzka et al. 2024, whose refit has a maximum amplitude of ~220 pc and a wavelength shrinking from ~4 kpc at Canis Major to ~0.4 kpc at Cygnus X).

**Data.** The Edenhofer et al. (2024) map on Zenodo (DOI 10.5281/zenodo.8187942, record 10658339): the mean map and its uncertainty, a 3.3 GB HEALPix FITS file, plus the paper's 2 kpc "validation" map (4.1 GB), needed for the far ends of the Radcliffe Wave (Canis Major, Cygnus X). The record holds 123 GB in all; its 12 posterior samples (19.5 GB) are one array, so single distance shells can be fetched by HTTP byte range (3.1 MB per sample and shell). The map is also available through the `dustmaps` package. Cloud distances and bridging sightlines (Zucker et al. 2020; Alves et al. 2020, Harvard Dataverse 07L7YZ, K16GQX); the Local Bubble model built from this map (O'Neill et al. 2024, Dataverse INB1RB) and Zucker et al. (2022)'s cluster tracebacks and tables (Dataverse E8PQOD, ZU97QD, 1VT8BC). Young star clusters from Hunt & Reffert (2023, CDS J/A+A/673/A114), since B8 is not built yet. Second map: Vergely et al. (2022), 10 pc cube on CDS (222 MB).

**Method + visualizations.**
- The keeper: a volume-rendered 3D fly-through (plotly or a WebGL viewer) starting at the Sun and passing Taurus and Orion.
- Face-on and edge-on maps of the dust density, with the walls of the Local Bubble.
- The Radcliffe Wave: fit a damped sine wave to the cloud positions and overlay it on the map.
- Young stars (younger than ~20 Myr, from B8) over the dust: star formation on the bubble's surface.
- Integrated extinction to any distance → a dereddening tool for B1 and B8.

**Check against.** Alves et al. 2020 (length 2.7 kpc, amplitude ~160 pc) and Konietzka et al. 2024 (Extended Data Table 2); Zucker et al. 2022 (the bubble's history); O'Neill et al. 2024 (the bubble's wall, from this same map: median distance 170 pc, thickness 35 pc, mass 6.0 × 10⁵ M☉); the paper's own Fig. 4, read back through its colour bar; a second map (Bayestar19 or Vergely et al. 2022).

**Watch out for.** Uncertain stellar distances smear dust along lines of sight. The map starts at 69 pc, so the very nearest clouds are missing. Maps from different teams disagree in detail; compare two. Rendering choices (opacity, colour) can create apparent structures. Interpolating the released *mean* map in log space (the authors' `interp2box.py`, `dustmaps`) reads 4-8% higher than averaging the interpolated posterior samples, which is what the paper's figures do (Hölder's inequality; measured on the samples). `dustmaps`' integrated mode places each voxel's cumulative extinction half a voxel too close, and the release's `interp2lbd.py` skips the radial interpolation (it uses the lower shell twice). Neighbouring voxels are positively correlated in the posterior (+0.3 to +0.5), so quadrature error sums are lower bounds. The map's ZGR23 extinction curve is closest to Cardelli R_V = 2.76, so A_V values from R_V = 3.1 work (e.g. Hunt & Reffert's clusters) run 11% higher for the same reddening. Alves et al. (2020)'s printed eq. (1) with their Extended Data Table 2 (δ = 6.68) never rises more than 61 pc, not 160 pc; use Konietzka et al.'s fully specified model.

**Builds on.** B1–B3 (Gaia), B8 (young clusters), H2 (gas), H5 (dust at all wavelengths), L3 (the supernovae that blew the bubble).

---

# I. Black holes

Black holes already run through the atlas: C1 hears two merge, C5 maps the merging population, H1 weighs the one at our Galaxy's centre. This section covers the rest, from ~10 M☉ black holes in X-ray binaries to billion-solar-mass ones powering quasars, to the first images of a black hole's shadow. Two numbers recur. The Schwarzschild radius is `R_s = 2GM/c² ≈ 3 km × (M/M☉)`. The Eddington luminosity, where radiation pressure balances gravity on infalling gas, is `L_Edd ≈ 1.26×10³⁸ (M/M☉) erg/s`.

**Suggested black hole arc:** I4 (stellar) → C5 (merging) → H1 (Sgr A*) → I5 + I6 (the shadow) → I1 (black holes and galaxies) → I2 + I3 (quasars).

## I1. The M–σ relation: black holes grow with their galaxies ★☆☆

**Science question.** A galaxy's central black hole weighs only ~0.5% as much as its central bulge, and its gravity dominates only a tiny region. Yet its mass tracks the random speeds σ of bulge stars hundreds to thousands of times farther out. How tight is that relation, and what does it imply?

**Background.** The M–σ relation (Ferrarese & Merritt 2000; Gebhardt et al. 2000) is roughly M ∝ σ^4–5 with ~0.3 dex scatter. The black hole's gravity matters only inside its "sphere of influence", `r ≈ GM/σ²`: tens of parsecs, against kiloparsec-sized bulges. The relation therefore implies that black holes and galaxies regulated each other's growth, probably through energy fed back from the black hole during its quasar phase.

**Data.** van den Bosch (2016) compilation of ~230 dynamically measured black hole masses with σ, luminosities and sizes (GitHub `remco-space/Black-Hole-Mass-compilation`, formerly `remcovandenbosch/`; VizieR `J/ApJ/831/134`). The Kormendy & Ho (2013) tables (ARA&A supplement; also in the arXiv source, Table 3 only as PostScript). The AGN Black Hole Mass Database (Bentz & Katz 2015) for masses from reverberation mapping (I2).

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

**Background.** The visible companion star's radial-velocity curve (period P, semi-amplitude K) gives the **mass function** `f(M) = P K³ / (2πG) = M_BH³ sin³i / (M_BH + M★)²`. This is a strict *lower limit* on the unseen mass. If f(M) exceeds ~3 M☉, above any possible neutron star, the object is a black hole whatever the geometry. V404 Cyg has f(M) = 6.08 ± 0.06 M☉. Cygnus X-1 (1972), the first black hole identified, was revised to 21.2 ± 2.2 M☉ in 2021 (Miller-Jones et al.). Known masses show a gap between ~2 and ~5 M☉ separating neutron stars from black holes. Gaia is now finding *dormant* black holes through the wobble of their companions, e.g. Gaia BH1, found in DR3's orbit solutions and weighed at 9.6 M☉ with follow-up velocities (DR3's orbit alone gives 12.8 ± 2.0).

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

**Background.** At low radio frequencies the old, extended lobes of radio galaxies are brightest. LOFAR's LoTSS DR2 imaged 27% of the northern sky at 144 MHz with 6″ resolution and catalogued 4,396,228 sources (Shimwell et al. 2022). A radio galaxy counts as a "giant" if it spans more than 0.7 Mpc. Machine learning plus citizen science raised the number of known giants above 11,000 (Mostert et al. 2024). Porphyrion (Oei et al. 2024, *Nature*) spans 7 Mpc, from a galaxy seen when the universe was ~6.3 Gyr old. The lengths are projections, so true sizes are larger. The population tells us how long jets stay switched on and how the cosmic web confines them.

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

# M. Binary and triple stars

About half of Sun-like stars have at least one companion, and most massive stars do. B9 used eclipsing pairs to weigh stars; C9, I4 and I11 found neutron stars and black holes in binaries; L5 chased the survivors of a pair that blew up. This section is about the companions themselves: how pairs and triples orbit, how often they form, how they trade mass and merge, and how the widest pairs test gravity itself.

Two stars obey Kepler's third law, `M₁ + M₂ = a³ / P²` (a in AU, P in years, M in solar masses), so any pair with a measured orbit and distance is a scale. Three bodies have no general solution. Stable triples are *hierarchical*: a close pair plus a distant third star, which slowly tugs the inner orbit around and can pump up its eccentricity (von Zeipel–Kozai–Lidov cycles).

**Suggested binary-star arc:** M1 (weigh stars by their orbits) → M2 (how common are companions?) → M4 (triples) → M5 (stars that touch and merge) → M6 (white dwarf pairs for LISA) → M3 (wide binaries as a test of gravity).

## M1. Weigh stars with binary orbits: Sirius B and its cousins ★☆☆

**Science question.** In 1844 Bessel noticed Sirius wobbling across the sky and deduced an unseen companion. That companion, Sirius B (first seen in 1862), turned out to be a white dwarf: the Sun's mass packed into a ball the size of the Earth. Can we fit visual binary orbits ourselves, weigh the stars with Kepler's third law, and watch the mass–luminosity relation emerge from hundreds of pairs?

**Background.** In a visual binary both stars are resolved, and the companion traces an ellipse around the primary: the projection of the true orbit onto the sky. Seven elements describe it (P, T, e, a″, i, Ω, ω). The angular semi-major axis a″ divided by the parallax ϖ gives a in AU, and then `M₁ + M₂ = (a″/ϖ)³ / P²`. Individual masses need each star's own wobble against the background (an astrometric orbit) or radial velocities (A3). The first binary orbit was computed by Savary in 1827 (published in the *Connaissance des Temps* for 1830), for ξ Ursae Majoris. Sirius A and B orbit every 50.13 years, with masses of 2.063 ± 0.023 and 1.018 ± 0.011 M☉ (Bond et al. 2017, from Hubble astrometry). α Centauri A and B, the nearest pair, orbit every 79.9 years.

**Data.** The Sixth Catalog of Orbits of Visual Binary Stars (ORB6, maintained with the Washington Double Star Catalog): a 1 MB text file of about 4,000 orbits graded from 1 (definitive) to 5 (indeterminate), at `https://www.astro.gsu.edu/wds/orb6/orb6orbits.txt` (https works again as of October 2026). Its documented column layout is not always respected: long periods in days and large errors spill into neighbouring columns, so parse by token overlap, not fixed columns (M1 found 34 periods misread). Gaia DR3 parallaxes and photometry, cross-matched by position; stars as bright as Sirius A are not in Gaia and need Hipparcos parallaxes. Sirius B (G = 8.5) *is* in DR3 (374.49 ± 0.23 mas), but its parallax is biased low by its orbital curvature during DR3; M1's scanning-law simulation corrects it to 379.0 ± 2.4 mas, matching Hipparcos (the test needs Gaia's 570 MB commanded scanning law). Bond et al. (2017) tables of Sirius B's positions.

**Method + visualizations.**
- Fit Sirius B's orbit from the published positions with our own least-squares (Thiele–Innes) fit plus MCMC → P, a″, e, i; with the Hipparcos parallax → the total mass, then the individual masses from the published astrometric orbit of Sirius A.
- An animated gallery of famous orbits from ORB6 elements, all to one scale: Sirius, Procyon, α Cen, ξ UMa, Castor, 61 Cyg.
- Dynamical masses for every grade 1–2 ORB6 orbit with a good Gaia parallax → total mass against total luminosity: the mass–luminosity relation from orbits alone, against B9's eclipsing-binary version.
- Sirius B as a white dwarf: its mass and radius (from temperature and luminosity) on B10's mass–radius diagram, with the Chandrasekhar curve. Its gravitational redshift is N3's showcase.

**Check against.** Bond et al. 2017: P = 50.13 yr, M_A = 2.063 M☉, M_B = 1.018 M☉. The mass–luminosity relation of Eker et al. (2018), as in B9.

**Watch out for.** Mass scales as (a″/ϖ)³, so a 3% error in the parallax or in a″ is a 9% error in mass. Many ORB6 orbits are grade 4–5, short arcs of long orbits, and useless for masses. Gaia DR3 parallaxes of close pairs can be biased by the companion's light; check RUWE. In practice only 55 of ORB6's 630 grade 1–2 orbits have a clean Gaia DR3 parallax (RUWE < 1.4, 5% precision): Hipparcos is still needed for most of them. When the observed arc is shorter than the period, P, a and e are degenerate.

**Builds on.** E2 and A3 (Keplerian orbits), B9 (masses from eclipses), B10 (white dwarfs), B2 (Gaia's nearby stars).

## M2. How common are companions? A census of binary stars ★★☆

**Science question.** Is the Sun, alone, typical? How often do stars have companions, how far apart are they, and does that depend on the star's mass? Can we rebuild the binary census from Gaia's 813,687 non-single stars and the classic spectroscopic catalogue, correcting for the companions each method cannot see?

**Background.** Among Sun-like stars within 25 pc, Raghavan et al. (2010) observed 56% single, 33% binary, 8% triple and 3% higher-order systems, and estimated that 54% are single after correcting for missed companions; for F and G dwarfs within 67 pc, Tokovinin (2014) finds 54:33:8:4:1 for systems of one to five or more stars. Their periods follow a broad log-normal peaking near 10⁵ days (tens of AU), from under a day to millions of years. The fraction with companions rises steeply with mass, from about a quarter of M dwarfs to most O stars (Moe & Di Stefano 2017). Close binaries have circular orbits because tides circularize them, up to a cut-off period that grows with age, from a few days in young clusters to 10–15 days in old ones (Meibom & Mathieu 2005): a clock in itself. There is an excess of "twins" with mass ratios above 0.95, out to separations of hundreds to thousands of AU (El-Badry et al. 2019), which hints that they grew together from a shared disk. Gaia DR3 flagged 813,687 non-single stars: astrometric, spectroscopic and eclipsing orbits, plus accelerating and wobbling sources.

**Data.** Gaia DR3 `nss_two_body_orbit` and the other non-single-star tables, joined to `gaia_source` for distances and colours (queried as in B1 and I11). SB9, the 9th Catalogue of Spectroscopic Binary Orbits (`sb9.astro.ulb.ac.be`, which timed out in October 2026; VizieR's copy B/sb9 has 4,079 systems and 5,099 orbits). Published tables of Raghavan et al. 2010 and Moe & Di Stefano 2017. DEBCat (B9).

**Method + visualizations.**
- Period–eccentricity diagrams for SB9 and Gaia's spectroscopic orbits: circular orbits at short periods, and the cut-off for young and old stars.
- Gaia's binaries on B1's HR diagram, including the twin sequence B1 found 0.75 mag above the main sequence.
- The period distribution within 100 pc, joining Gaia's spectroscopic and astrometric orbits (days to years) with resolved wide pairs (M3), each corrected for what its method can detect → compared with Raghavan's log-normal.
- Mass-ratio distributions and the twin excess.
- Multiplicity against primary mass, from M dwarfs to B stars.

**Check against.** Raghavan et al. 2010 (log-normal peak at log P = 5.03 days, σ = 2.28; 46% of Sun-like stars have companions); Moe & Di Stefano 2017; Gaia DR3 stellar multiplicity (Gaia Collaboration, Arenou et al. 2023).

**Watch out for.** Every method has a blind spot. Spectroscopy misses long periods and face-on orbits, Gaia DR3's orbits cover periods shorter than its 34-month span, and imaging misses close pairs. Correcting for completeness is the whole game. Gaia DR3 orbits passed strict filters, so they are a biased subset, and some are spurious (I11's warning). Unresolved companions hiding among the "single" stars push the single fraction up.

**Builds on.** B1 (the twin sequence), B9 (DEBCat), A3 (spectroscopic orbits), I11 (Gaia orbits and quality cuts).

## M3. Wide binaries: a test of gravity at low acceleration ★★★

**Science question.** Two Sun-like stars 10,000 AU apart pull on each other with only ~10⁻¹⁰ m/s². That is the regime where galaxy rotation curves seem to need dark matter (D4). MOND, the best-known alternative, says gravity itself strengthens there, so the widest pairs should orbit ~20% faster than Newton allows. Gaia has catalogued over a million wide binaries, and published analyses disagree about what they show. Can we see why?

**Background.** Below an acceleration a₀ ≈ 1.2 × 10⁻¹⁰ m/s², MOND (Milgrom 1983) boosts gravity; that reproduces the radial acceleration relation D4 measured (g† = 1.11 × 10⁻¹⁰ m/s²) without dark matter. Inside a binary dark matter is irrelevant, so wide binaries test the force law directly. We see only the sky-projected separation and velocity difference of each pair, so the test is statistical: compare the ratio of the projected velocity difference to the Newtonian circular speed at that separation, between tight pairs (high acceleration) and wide ones. Chae (2023, 2024) reported a boost of the size MOND predicts. Banik et al. (2024) analysed 8,611 pairs within 250 pc and preferred Newton at 19σ. Saglia et al. (2025) added precise radial velocities for a 3D test and found Newtonian motion; reanalyses have continued into 2026. Sample cuts and hidden third stars, which add velocity, are the main suspects.

**Data.** The El-Badry, Rix & Heintz (2021) catalogue of 1.3 million wide binaries within 1 kpc from Gaia EDR3, with chance-alignment probabilities: Zenodo record 4435257 (1.4 GB full table; a 354 MB "shifted" catalogue of fake pairs for calibrating chance alignments). Gaia DR3 radial velocities. Chae's samples and code are public (Zenodo 8416435, 10062232, 10987140; his base file carries Hwang et al. 2022's per-pair eccentricities, whose full catalogue is a Google Drive file linked from github.com/HC-Hwang). Banik et al.'s sample is available only on request, so it has to be rebuilt from Gaia DR3 with Pittordis & Sutherland's (2023) pair search (all 8.9 million sources with parallax > 10/3 mas and G < 20). Chae et al.'s (2026) 36 pairs with precise radial-velocity differences are tabulated in github.com/seratsaad/wb3d-gamma, with Saad & Ting's model code.

**Method + visualizations.**
- The catalogue on the sky and in distance; the separation distribution, and the chance alignments from the shifted catalogue.
- Projected velocity difference against separation: the Keplerian v ∝ s^(−1/2) envelope for each total mass, with MOND's prediction, including the Galaxy's external field, drawn on.
- A forward model: random orientations and orbital phases, an eccentricity distribution, hidden companions and line-of-sight interlopers, calibrated on tight pairs, then predicted for wide ones under Newton and under MOND.
- The test itself, then repeated with each published sample's cuts. Which choices move the answer?

**Check against.** Banik et al. 2024: α_grav = −0.021 (+0.065/−0.045), where 0 is Newton and 1 is MOND (MOND excluded at 16σ); Chae 2024 for the opposite conclusion; El-Badry et al. 2021's counts.

**Watch out for.** Undetected inner binaries inflate velocity differences and mimic a boost. Their rate is the dominant uncertainty, so fit it rather than assume it. Projection and eccentricity matter as much as the force law. Gaia's proper-motion errors must be propagated pair by pair, and perspective effects removed for nearby pairs. El-Badry et al.'s search keeps only pairs slower than a 5 M☉ circular orbit (+2σ), truncating ṽ near √(5 M☉/M), and its R < 0.01 filter removes fast wide pairs too; Chae's 26,615 needs two cuts found only in his code (dec > −28° at all distances, a colour–magnitude line). In the (r_sky, ṽ) forward model the eccentricity law trades off against gravity (left free, γ ≈ 3 buys α ≈ 0.1), so constrain it with the pairs' own v–r angles, modelling their companions. Test any method on mock Newtonian and MOND skies first: Chae's statistic reads ~1.5 on a Newtonian mock of his own sample. This reproduces a live debate; do not present it as a verdict.

**Builds on.** D4 (the radial acceleration relation), M2 (multiplicity), M4 (triples), B2 (Gaia astrometry).

## M4. Three suns: triple stars and the three-body problem ★★☆

**Science question.** Tatooine had two suns; many real stars have three. How do triples stay stable when the three-body problem has no general solution? And how does a distant third star reshape a close pair, sometimes until it stops eclipsing?

**Background.** Stable triples are hierarchical: an inner binary with a third star orbiting the pair from much farther out. Mardling & Aarseth (2001) give the stability limit: the outer pericentre must be roughly three or more times the inner orbit's size, depending on the masses and eccentricities. The third star makes the inner orbit precess and, if the two orbits are tilted by more than ~39°, trades inclination for eccentricity in slow von Zeipel–Kozai–Lidov cycles. Combined with tides, these cycles shrink binaries: 96% of spectroscopic binaries with periods under 3 days have a third companion (Tokovinin et al. 2006; Fabrycky & Tremaine 2007). Precession also tilts eclipsing pairs: SS Lacertae eclipsed early in the 20th century and had stopped by about 1950 (Torres & Stefanik 2000). Showcases: α Centauri, where Proxima orbits the A–B pair every ~550,000 years at ~8,700 AU (Kervella et al. 2017); HD 181068, a "trinary" seen by Kepler (Derekas et al. 2011); and TIC 290061484, three stars of 6–8 M☉ that all eclipse one another, with an outer orbit of only 24.5 days around a 1.8-day pair (Kostov et al. 2024). The triple pulsar PSR J0337+1715 tests gravity itself (N7).

**Data.** Tokovinin's Multiple Star Catalog (MSC): about 2,000 hierarchies of 3–7 stars, on VizieR (`J/ApJS/235/6`) and updated online at `www.ctio.noirlab.edu/~atokovin/stars/` (plain http). TESS light curves of TIC 290061484 and Kepler light curves of HD 181068 and KOI-126 via `lightkurve` (A1). Published orbits of α Cen AB and Proxima (Kervella et al. 2016, 2017).

**Method + visualizations.**
- The MSC's architectures: outer against inner period, with the Mardling–Aarseth stability line, and mutual inclinations where known.
- Eclipse timing variations: eclipse times of TIC 290061484's inner pair measured from TESS → the light-travel-time and dynamical delays that reveal the 24.5-day outer orbit.
- REBOUND integrations of real triples from their published elements (labelled simulations): Kozai–Lidov cycles, the precession that switched off SS Lac's eclipses, and what happens when a triple is pushed past the stability limit (chaos, then ejection).
- Proxima's orbit around α Cen AB, drawn to scale with the Solar System and the Oort cloud.

**Check against.** Kostov et al. 2024 (P_out = 24.5 d, P_in = 1.8 d, masses 6.85, 6.11 and 7.90 M☉); Kervella et al. 2017 (Proxima: P ≈ 547,000 yr, e ≈ 0.50).

**Watch out for.** Eclipse times need careful templates, because spots and blended light shift them. Kozai cycles need general-relativistic and tidal precession switched on, since either can quench them (REBOUNDx provides both). An integration of a chaotic triple is not a prediction for that system, only a statistic. The MSC is a compilation with uneven completeness.

**Builds on.** M1 (visual orbits), B9 (eclipsing binaries), A1 (light curves), E1 (REBOUND), I11 (hidden triples that mimic dark companions), M3 (triples are the wide-binary systematic).

## M5. Stars that touch, trade mass and merge ★★☆

**Science question.** In a close binary the stars can swell until they touch. Then gas flows from one to the other, and sometimes the pair spirals into a single star. In 2008 astronomers watched one do it: V1309 Scorpii, a contact binary whose period had been shrinking for years, erupted as a "red nova". Can we find contact binaries by the hundred thousand, follow that merger, and learn from a famous merger prediction that failed?

**Background.** Each star in a binary sits in a teardrop-shaped region, its Roche lobe, inside which its own gravity dominates. A star that outgrows it spills gas onto its companion through the inner Lagrange point L1. In **Algol** systems the less massive star is the more evolved one, a puzzle (the "Algol paradox") until astronomers realised it had given most of its mass away. **Contact** (W Ursae Majoris) binaries share one envelope; they vary continuously, with periods mostly under a day and a sharp lower cut-off near 0.22 days (Rucinski 1992). Their periods correlate with luminosity tightly enough to measure distances. Losing angular momentum shrinks them until they merge. V1309 Sco's ~1.4-day period fell by about 1% between 2001 and 2007, its light curve changed shape, and it erupted in 2008 (Tylenda et al. 2011). Molnar et al. (2017) predicted that KIC 9832227 would merge around 2022; Socia et al. (2018) traced its apparent period decrease to one erroneous early timing, and the star is still there. T Coronae Borealis, a white dwarf feeding on a red giant, erupts as a nova roughly every 80 years (1866, 1946); it had not erupted by September 2026.

**Data.** The ZTF catalogue of 781,602 periodic variables, with contact (EW) binaries among its classes (Chen, Wang & Deng 2020; VizieR), and Gaia DR3's eclipsing-binary candidates (`vari_eclipsing_binary`, as in B9; also on VizieR as I/358/veb; note that Gaia's period search stopped at 0.2 d). Kepler light curves of KIC 9832227 via `lightkurve`, plus TESS (2019–2024) and the older survey photometry needed to check the merger prediction: NSVS 1999–2000 (VizieR II/287), SuperWASP 2004–2010 (NASA Exoplanet Archive DR1; VeSPA), ASAS-SN 2014–2018 and AAVSO 2015–2026. V1309 Sco: the OGLE-III/IV photometry 2001–2010 published with Tylenda et al. (2011) is public at `astrouw.edu.pl/ogle/ogle3/V1309_SCO/` (*corrected 2026-10-09: an earlier version of this spec said it was not in a public archive*). T CrB: the AAVSO International Database, readable without an account through the VSX API (`aavso.org/vsx/index.php?view=api.object&ident=T+CrB&data=N&fromjd=..&tojd=..&format=json`; *corrected 2026-10-09: an earlier version said downloads need an account*); the server stops at 30 s per request, so 2024–26 needs windows of a week or two.

**Method + visualizations.**
- Roche geometry: equipotential surfaces and L1 for mass ratios of 1 and 0.1, in 3D.
- Contact binaries in ZTF and Gaia: the period distribution and its cut-off; the period–luminosity relation with Gaia distances → a distance indicator.
- Eclipse timing (O − C) for KIC 9832227 from Kepler: the period change that was not there, and how one wrong time faked it.
- V1309 Sco from Tylenda et al.'s values: the shrinking period and the changing light curve before the merger.
- T CrB: its 1946 outburst and its recent pre-eruption behaviour in AAVSO data; if it erupts, its decline.

**Check against.** Tylenda et al. 2011; Socia et al. 2018; the contact-binary period–luminosity relations of Chen et al. (2018).

**Watch out for.** Contact, semi-detached and ellipsoidal variables look alike in survey photometry. O − C diagrams mix real period changes with light-travel time from third stars (M4) and with spots. Period–luminosity relations differ for early- and late-type contact binaries. Merger predictions from period changes are fragile: that is the KIC 9832227 lesson. Time conventions: the NSVS archive labels its time column "JD − 2400000" but stores MJD (JD − 2400000.5), and "MJD" means JD − 2400000.0 in the Kepler EB catalogue; that half day plus one miscounted orbit is what produced the 2022 prediction. AAVSO's unvalidated visual reports include false "eruptions" of T CrB (several at magnitude 4–5 in 2025–26); compare with same-night CCD data.

**Builds on.** B9 (eclipsing binaries), A1 (Kepler light curves), M4 (third bodies in O − C), I4 and L5 (mass transfer onto compact stars), B10 (white dwarfs).

## M6. White dwarf pairs: LISA's guaranteed sources ★★☆

**Science question.** ZTF J1539+5027 is a pair of white dwarfs that orbit each other every 6.91 minutes and eclipse on every orbit. Gravitational waves are draining their orbital energy, so the eclipses arrive measurably earlier each year. Can we measure that from public ZTF data, and predict what LISA, the space-based gravitational-wave observatory, will hear from such pairs?

**Background.** A binary of orbital period P emits gravitational waves at a frequency of 2/P. For orbits of minutes that is millihertz: too low for LIGO, right for LISA (adopted by ESA in January 2024, launch planned for about 2035). The orbital decay follows C9's quadrupole formula, so Ṗ depends only on the chirp mass and P. For ZTF J1539+5027, Ṗ = (−2.373 ± 0.005) × 10⁻¹¹ s/s (Burdge et al. 2019); Burdge's updated ephemeris, from six years of high-speed photometry, gives (−2.36812 ± 0.00044) × 10⁻¹¹ (published in Shaw et al. 2024, AJ). Kupfer et al. (2024) analysed 55 candidates: 40 should be detected within LISA's four years, 18 of them ("verification binaries") within the first three months, with signals predictable in advance from telescope data. The Galaxy's many millions of unresolved white dwarf pairs will add up to a foreground hum.

**Data.** ZTF difference-image cutouts and their PSFs from IRSA's image archive, for forced photometry: IRSA's light-curve service lists only detections, which drop the faint phases of a star at the detection limit. That is ~5,900 small files; IRSA stages ~18% of them from slow storage (3-20 minutes each), and 5% of exposures have no difference image. Public PTF Level 1 images (2009-2012) extend the baseline to 17 years. Burdge et al. 2019 did not tabulate their eclipse times; their ephemeris and its update are the check. The Kupfer et al. (2024) list of LISA binaries with their Gaia DR3 distances (Tables 1-2, parsed from the arXiv LaTeX source). The analytic LISA sensitivity curve of Robson, Cornish & Liu (2019).

**Method + visualizations.**
- ZTF J1539's survey photometry folded on 6.91 minutes; fit P and Ṗ from the drift of the eclipse phase across seasons (an O − C parabola like C9's Hulse–Taylor plot) → the chirp mass, compared with Burdge et al. and with GR's prediction from their component masses.
- The quadrupole formula for every Kupfer et al. binary → frequency and strain; the LISA sensitivity curve with the verification binaries on it, and their signal-to-noise after 4 years.
- Time to merger for each pair, and their likely fates: AM CVn mass transfer, a merger into one star, or a Type Ia supernova (L5).
- The Galactic foreground: a simple population model scaled from the white-dwarf density of the 40 pc sample (O'Brien et al. 2024; B10 not yet built), with Maoz et al.'s statistics of white-dwarf pairs, as in Korol et al. (2022).

**Check against.** Burdge et al. 2019: P = 6.91 min, Ṗ = −2.373 × 10⁻¹¹ s/s; Shaw et al. 2024 (Burdge's update): Ṗ = −2.36812(44) × 10⁻¹¹; Kupfer et al. 2024: 18 verification and 40 detectable binaries; Korol et al. 2022: ~60,000 resolvable pairs for 2.7 × 10⁹ disc white dwarfs.

**Watch out for.** The target is faint for ZTF and its 30 s exposures smear a 6.9-minute orbit, so model the exposure time; timing is coarse per exposure, but over the full archive the eclipse time comes out good to ~0.3 s and Ṗ to 0.2%. Burdge et al.'s high-speed photometry is the precise check. Barycentric time corrections are essential, and the ZTF headers have two traps: `OBSJD` is ~0.7 s before the shutter opens, and "HJD" is not heliocentric (use `SHUTOPEN`/`SHUTCLSD`). With free eclipse depths a light-curve fit can lock on half an orbit away. LISA SNRs at 1-3 mHz in the first months depend mostly on the assumed Galactic foreground. Tides and mass transfer also change periods, so for the tightest pairs Ṗ is not pure gravitational radiation. LISA's sensitivity depends on mission assumptions; use the published analytic curve and say so.

**Builds on.** C9 (orbital decay), C1 (gravitational waves), B10 (white dwarfs), L5 (white dwarf mergers), A1 (folding light curves).

---

# N. General relativity

Newton's gravity is a force. Einstein's (1915) is geometry: mass and energy curve spacetime, and freely falling bodies and light follow the straightest possible paths through it (geodesics). The curvature is encoded in the **metric**, which gives the interval between nearby events. Schwarzschild found the metric outside a non-spinning mass in 1916, Kerr the spinning one in 1963. GR has already appeared throughout this atlas: C1 and I7 heard spacetime ripple, C9 watched a binary decay, H1 measured S2's gravitational redshift, I5 and I6 saw a black hole's shadow. This section tests the theory where it started, with planets, starlight, clocks and gyroscopes, and works with its mathematics directly: the symmetries (Killing vectors) that make orbits solvable, and a fall into a black hole.

In weak gravity, departures from GR are summarised by the **parametrized post-Newtonian (PPN)** parameters: γ measures how much space curvature a unit mass produces (light bending and time delay measure it), and β how nonlinear gravity is (perihelion precession). GR has γ = β = 1.

**Suggested relativity arc:** N7 (the scorecard) → N2 (light bends) → N3 (clocks) → N1 (orbits precess) → N4 (gyroscopes) → N6 (symmetries) → N5 (falling in).

## N1. Precessing orbits: from Mercury to S2 ★★☆

**Science question.** Mercury's orbit turns slowly: its perihelion advances by ~575″ per century. Newtonian tugs from the other planets explain all but 43″, and in 1915 Einstein explained those 43″, the first success of GR. Can we reproduce that budget with an N-body integration from real ephemerides, then follow the same effect up a ladder of ever stronger gravity: binary stars, binary pulsars and the star S2 at the Galactic Centre?

**Background.** In GR an elliptical orbit does not close. Each orbit advances by `Δω = 6πGM / (c² a (1 − e²))`: 42.98″ per century for Mercury, 8.6″ for Venus, 3.8″ for Earth. Le Verrier found the anomaly in 1859 (38″; Newcomb later 43″); Vulcan, the planet proposed to explain it, does not exist. MESSENGER's radio ranging separated GR's share from the Sun's oblateness (J2☉), giving β − 1 = (−2.7 ± 3.9) × 10⁻⁵ (Park et al. 2017). Up the ladder: eccentric eclipsing binaries precess through tides, rotation and GR together. DI Herculis turned at only about a quarter of the predicted rate for decades, until its stars' spin axes were found tipped almost sideways (Albrecht et al. 2009). PSR B1913+16's periastron advances by 4.23° per year (C9). S2's orbit turns by 12′ per orbit, which GRAVITY measured in 2020 (f_SP = 1.10 ± 0.19, where GR is 1); H1 found that older speckle positions fake a much larger value.

**Data.** JPL Horizons state vectors and osculating elements for the planets (E2's queries; DE440/DE441); N1 used DE440's own SPICE file (de440.bsp, 120 MB, 1550–2650) for initial conditions, masses and the answer key, and Horizons' DE441 for 1000 BC–AD 5000. REBOUND plus REBOUNDx 5.1, which supports REBOUND 5 and includes `gr`, `gr_potential` and `gr_full` effects (source-only on PyPI; it compiles on install with the system gcc). TESS apsidal-motion measurements for eccentric eclipsing binaries (Baroch et al. 2021, "a test of gravitational theories") and DI Her's published parameters. GRAVITY 2020 and H1's S2 fits; C9's pulsar timing parameters.

**Method + visualizations.**
- Integrate the Solar System for a few centuries from Horizons initial conditions, first Newtonian only, then with REBOUNDx's GR term. Mercury's perihelion longitude against time in each case → ~532″/cy from the planets plus 43″/cy from GR, and each planet's share as a bar chart (Venus, Jupiter and Earth dominate).
- Check against Horizons' own elements over the same centuries (DE441 includes GR), and against Le Verrier's and Newcomb's historical numbers.
- An exaggerated animation: Mercury with GR turned up a million times, tracing a rosette.
- The precession ladder: Mercury, Earth, the TESS binaries, DI Her, the Hulse–Taylor and double pulsars, and S2 on one plot of precession per orbit against GM/(c²a), with GR's line through them.

**Check against.** 42.98″/cy from GR and ~532″/cy from the planets; Park et al. 2017; GRAVITY Collaboration 2020 (f_SP = 1.10 ± 0.19); Albrecht et al. 2009.

**Watch out for.** Precession depends on the reference frame. In the Earth's precessing equatorial frame Mercury's perihelion moves ~5,600″/cy, mostly the precession of the equinoxes. Even fixed frames disagree: the textbook 575″, ~532″ and 43″/cy are rates *within Mercury's own orbital plane* (Park et al. 2017 use Mercury's mean orbit plane); the longitude of perihelion Ω + ω in the fixed J2000 ecliptic advances only 571.9″/cy (planets 528.9″/cy) and on the ICRF equator 560.9″/cy, because the orbit plane slides along the reference plane (N1 found this). The rate at J2000 also depends on the fitting window at ±0.02″/cy: Park et al. fitted AD 1000–3000 with a quadratic plus 14 periodic terms. Keep numerical drift far below 1″/cy (IAS15, or a small time step) and check by integrating backward. The Sun's oblateness adds ~0.03″/cy. In stars, apsidal motion is dominated by tides and rotation, and GR's part comes out only with good stellar-structure models.

**Builds on.** E2 (Horizons), E1 (REBOUND), H1 (S2), C9 (binary pulsars), B9 (eclipsing binaries).

## N2. Bending starlight: from the 1919 eclipse to Cassini ★☆☆

**Science question.** On 29 May 1919 two British expeditions photographed the stars around the eclipsed Sun and found them pushed outward by about the 1.75″ Einstein predicted, twice the Newtonian value. The result made Einstein famous overnight, and it has been argued over ever since. Can we redo the measurement from the expeditions' own published star positions, then follow light bending and its partner, the Shapiro delay, to today's precision of 10⁻⁵?

**Background.** Light passing the Sun at impact parameter b is deflected by `α = (1 + γ)/2 × 4GM/(c²b)`, which is 1.75″ × (R☉/b) in GR (γ = 1). Treating light as Newtonian particles gives half that (Soldner 1801: 0.84″). Stars near the eclipsed Sun therefore appear pushed radially outward. Between the eclipse plates and comparison plates taken months later at night, the images also shift, rotate and change scale, so these plate constants must be fitted together with the deflection. The Sobral 4-inch lens (7 plates, 7 stars) gave 1.98 ± 0.12″ and Príncipe 1.61 ± 0.30″ (probable errors; one sigma 0.18″ and 0.44″). The Sobral astrographic plates gave 0.93″, or 1.52″ with the scale change taken from refraction alone; the Sun's heat had warped the astrographic's coelostat mirror, and the report set its result aside (Dyson, Eddington & Davidson 1920). Critics later called that choice biased; a 1978 remeasurement of the surviving Sobral plates (published 1979) found 1.90 ± 0.11″ for the 4-inch and 1.55 ± 0.34″ for the astrographic, both consistent with Einstein (Harvey 1979; Kennefick 2009). The 1922 Wallal eclipse, with many more stars, gave 1.72 ± 0.11″ (probable error) from its 15-foot cameras (Campbell & Trumpler 1923), and 1.75 ± 0.09″ with its five-foot cameras added (1928). Radio waves bend too: geodetic VLBI of thousands of radio sources all over the sky (1979–2010) gives γ = 0.99992 ± 0.00012 (formal error; Lambert & Le Poncin-Lafitte 2011). The same curvature delays light: the frequency shift of radio signals to the Cassini spacecraft, whose path passed 1.7 R☉ from the Sun's centre in June 2002, gave γ − 1 = (2.1 ± 2.3) × 10⁻⁵ (1σ; Bertotti, Iess & Tortora 2003). Cassini measured the delay's rate of change, not the delay itself. Gaia corrects every star position for the Sun's bending, which is still 4 mas at 90° from the Sun.

**Data.** Dyson, Eddington & Davidson 1920 (*Phil. Trans. R. Soc. A* 220, 291; free online), with tables of measured star displacements for each plate. The report prints full per-plate measures only for the 4-inch and Príncipe plates; for the astrographic, only per-plate summaries and the declination shifts of 5 stars. Campbell & Trumpler's Wallal tables (*Lick Observatory Bulletin*, on ADS). Cassini's raw tracking data: PDS Atmospheres Node, volumes CORS_0021–0028 (closed-loop Doppler, ~1.3 GB). Published VLBI, Viking, Cassini and Gaia values for the history plot.

**Method + visualizations.**
- Refit the 1919 plates: for each plate, the measured displacements against star positions, with a linear plate model (the report's six constants: shift, scale and rotation in each coordinate) plus a deflection ∝ 1/r → α at the solar limb with its error. Repeat for Wallal.
- The famous diagram redrawn: each star's displacement as an enlarged arrow around the Sun's disk, with GR's and the half-deflection predictions.
- Deflection against distance from the Sun (1/r), with both expeditions' stars.
- A century of γ: |γ − 1| from 1919 to Cassini on a log scale, as in Will's review.
- The Shapiro delay: the extra round-trip time to Cassini against impact parameter during the 2002 conjunction (about a quarter of a millisecond near the Sun), and the frequency shift Cassini actually measured.

**Check against.** Dyson et al. 1920 (1.98″ and 1.61″); Harvey 1979 (1.90 ± 0.11″ and 1.55 ± 0.34″); Campbell & Trumpler 1923 (1.72 ± 0.11″) and 1928 (1.75 ± 0.09″); Bertotti et al. 2003; Anderson, Lau & Giampieri 2004 (Cassini's closed-loop data).

**Watch out for.** Plate constants and deflection are correlated, so the error depends on how many stars span the field; Sobral had 7. Changes in the telescope between the eclipse and comparison exposures (the Sun's heat warped the astrographic's coelostat mirror and blurred its images, while its focusing scale stayed unchanged; Dyson et al. 1920, p. 309) are systematic, not random. Report a fit with all plates and with the 1920 selection, and let the reader judge.

**Builds on.** E8 (eclipses), I6 and J5 (light bent by black holes and galaxies), C9 (Shapiro delay in pulsars), B2 (Gaia astrometry).

## N3. Clocks run slower downhill: gravitational redshift from GPS to white dwarfs ★★☆

**Science question.** GR says clocks deeper in gravity tick slower. GPS satellite clocks gain 38 microseconds a day on clocks on the ground, and positions would drift by ~11 km a day if engineers ignored it. In 2014 two Galileo navigation satellites were launched by mistake into eccentric orbits; their hydrogen-maser clocks now speed up and slow down by ~0.4 µs on every orbit as they climb and fall. Can we see that in public navigation data, and then measure the same effect on the Sun and on white dwarfs?

**Background.** A clock at gravitational potential Φ moving at speed v runs at a rate `1 + Φ/c² − v²/2c²` relative to a distant clock. For GPS that is +45.7 µs/day (higher up) and −7.2 µs/day (orbital speed); the two cancel for a satellite about 3,200 km up. In an eccentric orbit the total oscillates once per orbit, by `Δt = −(2/c²) √(GMa) e sin E` (E is the eccentric anomaly), about ±0.37 µs for Galileo GSAT-0201 and GSAT-0202. Delva et al. (2018) and Herrmann et al. (2018) fitted that modulation in about three years of clock data: (+0.19 ± 2.48) × 10⁻⁵ of GR's prediction, 5.6 times better than Gravity Probe A's 1976 rocket. Pound and Rebka (1960) first measured the effect over the 22.5 m height of a Harvard tower. On the Sun the redshift is 633 m/s, measured to about 1% with HARPS and a laser frequency comb in sunlight reflected by the Moon (González Hernández et al. 2020). On a white dwarf it is tens of km/s: Sirius B's lines are shifted by 80.65 ± 0.77 km/s (Joyce et al. 2018).

**Data.** Precise multi-GNSS orbits (SP3) and satellite clocks: ESA's navigation office products (open access at `navigation-office.esa.int/products/gnss-products/{GPS week}/`), for Galileo E18 (GSAT-0201), E14 (GSAT-0202) and normal Galileo satellites as controls. The 5-min clocks inside the SP3 files equal the 30-s clock files at those epochs, so the ~1 MB/day SP3 files suffice (the clock files are ~8 MB/day). ESA's series (ESA0MGNFIN, a 2023 reprocessing to 2022 and operational after, then ESA0OPSFIN from mid-2026) includes E14 and E18 only from August 2016; ten years stream ~3.4 GB, ~0.45 GB kept. ESA's weekly SLR residual summaries (to January 2025) check the orbits. Sirius B: Joyce et al. 2018 and M1's mass. Published Pound–Rebka, Gravity Probe A and GREAT results.

**Method + visualizations.**
- The GPS budget: the gravitational and velocity terms against orbital radius, crossing at ~1.5 Earth radii.
- The Galileo clocks: compute the periodic relativistic term from the SP3 orbits. IGS-convention clock products have it removed, so add it back to see the ~0.4 µs swing once per orbit. Then fit a scale factor on that term, plus a polynomial for clock drift → the fractional deviation from GR, with an honest error from a year of data.
- The eccentric satellites beside a normal one: the swing appears only in the eccentric ones.
- Sirius B: the predicted redshift GM/(Rc) from M1's mass and its radius, against the measured 80.65 km/s.
- The redshift ladder: tower, rocket, satellites, the Sun, white dwarfs, S2 (H1) and a neutron star's surface (z ≈ 0.2–0.3, C8), as fractional shift against GM/(rc²).

**Check against.** Delva et al. 2018: (+0.19 ± 2.48) × 10⁻⁵; Herrmann et al. 2018; Joyce et al. 2018: 80.65 ± 0.77 km/s; GPS: 38 µs/day (Ashby 2003, *Living Reviews in Relativity*).

**Watch out for.** Know the convention: IGS satellite clocks leave out the periodic relativistic term that users are meant to add, so the products already assume GR, and the test is whether any residual remains. Clocks are estimated together with orbits, and orbit errors at the orbital period (solar radiation pressure, eclipse seasons) leak into the clocks at exactly the signal's period; that was the GREAT teams' main systematic. Maser drifts and frequency jumps must be modelled. A white dwarf's redshift is mixed up with its own radial velocity; Sirius B's orbit supplies that. Found when this was built: E18's clock carries a ~160 ps once-per-orbit signal locked to the Sun (thermal or radiation pressure), which fakes α of ~10⁻³ in any single year, so Sun-locked terms must be fitted and the error calibrated with null templates; with them the test is systematics-limited at ~4 × 10⁻⁵, no tighter than the published ones despite three times more data. ESA's SLR pass residuals disagree from station to station at the millimetre level, so they bound the orbit term rather than correct it. The IGS convention also leaves J2 and tidal terms (110–150 ps, twice per orbit) in the clocks. E14's maser degrades from 2024.

**Builds on.** M1 (Sirius B), B10 (white dwarfs), H1 (S2's redshift), F4 (the solar spectrum), C8 (neutron stars).

## N4. Gyroscopes in curved space: geodetic precession and frame dragging ★★☆

**Science question.** A perfect gyroscope orbiting the Earth slowly turns for two reasons: because space around the Earth is curved (geodetic precession, 6.6″ per year in Gravity Probe B's orbit), and because the spinning Earth drags spacetime around with it (frame dragging, 0.039″ per year). Gravity Probe B measured both. In July 2026 laser-ranged satellites measured frame dragging to ~0.2%. Can we compute these effects ourselves, compare them with every measurement, measure frame dragging ourselves from public laser-ranging orbits, and find the same physics in binary pulsars whose spin axes wobble?

**Background.** The geodetic effect (de Sitter 1916) turns a gyroscope's axis in the orbital plane at `Ω = (3/2) (GM)^(3/2) / (c² r^(5/2))` for a circular orbit. Frame dragging (Lense & Thirring 1918) depends on the central body's angular momentum J. Gravity Probe B (2004–2005) carried four quartz spheres, the roundest objects ever made, in a 642 km polar orbit. It measured a geodetic drift of −6,601.8 ± 18.3 mas/yr against GR's −6,606.1, and frame dragging of −37.2 ± 7.2 mas/yr against −39.2 (Everitt et al. 2011); the Sun's own geodetic effect (+7.3 mas/yr north-south, −16.2 west-east at the guide star) and the star's proper motion were removed first. Frame dragging also turns orbits: a satellite's orbital plane (its node) drifts by ~31 mas/yr for LAGEOS, while the Earth's equatorial bulge makes it drift by over 100° per year, ten million times more. LARES 2 (launched 2022) was placed in an orbit with LAGEOS's supplementary inclination, so the two bulge-driven drifts cancel in the sum while frame dragging adds; Ciufolini et al. (2026, *Nature*) found μ = 1.001 ± 0.001 (stat) ± 0.002 (syst) of GR's prediction. Pulsars are gyroscopes too: pulsar B of the double pulsar precesses at 4.77 (+0.66/−0.65)° per year (Breton et al. 2008), refined to 5.16 ± 0.33 (Lower et al. 2024), against GR's 5.07; PSR J1906+0746 at 2.17 ± 0.11 against 2.23 (Desvignes et al. 2019); PSR B1913+16's pulse profile changes as its beam precesses (Kramer 1998), though its rate has not been measured independently; and PSR J1141−6545's orbit is tilted by its fast-spinning white dwarf companion (Venkatraman Krishnan et al. 2020).

**Data.** The ILRS analysis centres' weekly precise orbits of LAGEOS and LARES 2 (SP3, centimetre-level, Earth-fixed SLRF2020, 2-min sampling), eight centres, at the EUROLAS Data Center by anonymous FTP (`ftp://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/{sat}/{year}/{week}/`; CDDIS needs an Earthdata login). IERS EOP 20 C04 from the Paris Observatory (`hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now`; the IERS data centre's "latestVersion" file stopped in May 2025). GOCO06s from ICGEM, FES2004 ocean tides from the IERS Conventions, UT/CSR's monthly C20 for comparison. Published results: Everitt et al. 2011 and 2015 (GP-B), Ciufolini et al. 2026 (LAGEOS + LARES 2), the pulsar papers above, lunar laser ranging (Hofmann & Müller 2018). Two-line elements (CelesTrak; Space-Track's archive needs a free account) only for comparison.

**Method + visualizations.**
- GP-B's predicted geodetic and frame-dragging rates, integrated along its orbit with the Earth's J2, the Sun and the Moon (not just the circular formula) → −6,597 to −6,620 mas/yr depending on the orbit's mean radius, and −39.1.
- The node drifts of LAGEOS and LARES 2 from the weekly orbits: +125.1 and −125.3° per year, their sum −0.19° per year.
- Frame dragging measured: refit every weekly orbit with an own IERS-2010 force model (frame dragging off), measure the weekly "jumps" of each orbital plane beyond the model, fit their sum for μ with periodic terms where tide- and radiation-pressure-model errors appear, and the Earth's J2 from their difference.
- A diagram of the two effects on a gyroscope: the geodetic turn in the orbital plane, frame dragging across it.
- The spin-precession ladder: the Moon's orbit, GP-B, PSR J0737−3039B, J1906+0746, B1534+12, and the frame-dragging measurements, measured against predicted.

**Check against.** GP-B (Everitt et al. 2011): −6,601.8 ± 18.3 and −37.2 ± 7.2 mas/yr against −6,606.1 and −39.2; LARES 2 + LAGEOS: μ = 1.001 ± 0.001 ± 0.002 (Ciufolini et al. 2026); Breton et al. 2008 and Lower et al. 2024; CSR's C20 series.

**Watch out for.** We are not redoing GP-B's analysis: unexpected torques from electrostatic patches on its gyroscopes took five years to model, and its raw data (~1 TB at NSSDCA by request) are not easily usable. GP-B's geodetic prediction scales as r^−5/2, 2.3 mas/yr per km of mean orbital radius, and the published orbit numbers disagree by 9 km. Two-line elements are mean elements with kilometre-level errors (265–340 m cross-track): fine for the bulge-driven drift, useless for frame dragging. Do not fit the radiation-pressure scale C_R per weekly arc: in eclipse seasons it absorbs 0–19% of frame dragging (synthetic test) and fakes inclination drifts (−21 mas/yr in both satellites); hold it at one value per satellite. The de Sitter precession (~17 mas/yr per node) moves both planes like frame dragging and cannot be separated in this geometry; keep it at GR's value. Fit periodic terms at K1, K2, S2, P1, the radiation-pressure argument (node − Sun) and the year; K1's ~1,050-day period is partly degenerate with a trend over a few years. The IERS FES2004 file's 18.6-year equilibrium ocean tide comes out in anti-phase with the solid tide; with it, J2 from the pair is offset from CSR's series by 1.1 × 10⁻¹⁰ in C̄20 (frame dragging is unaffected). One EDC file (jcet LAGEOS 2025-04-05) is a 98-byte header; EDC refuses more than ~4 simultaneous FTP connections. The frame-dragging error budget is dominated by tide-model errors (K1) and by which analysis centre's orbits are used.

**Builds on.** N1 (orbital precession), C9 (the double pulsar), C2 (ATNF catalogue), N3 (GNSS orbits and frames), N5 and N6 (frame dragging near spinning black holes), E2 (orbits).

## N5. Falling into a black hole (labelled simulation) ★★☆

**Science question.** What would you see, and feel, falling into a black hole? Would you be stretched into spaghetti before reaching the horizon, or cross it without noticing? How long would the trip take by your own watch, and how would it look to a friend watching from far away? Can we compute it from GR for the real black holes this atlas has weighed, with the real sky behind them?

**Background.** Outside a non-spinning black hole spacetime is described by the Schwarzschild metric. Its time coordinate is a distant observer's time, in which the faller never quite reaches the horizon: their light reddens and fades, the redshift growing by a factor e every 4GM/c³ (85 s for Sgr A*). By the faller's own clock the horizon is reached in finite time, and locally nothing special happens there. Inside, the singularity is not a place but a moment in the future, reached within at most πGM/c³ of proper time: 66 s for Sgr A*, about a day for M87*. Tidal stretching scales as M/r³, so at the horizon it goes as 1/M². A 10 M☉ black hole stretches a person with ~10⁷ g long before they reach the horizon, while at Sgr A*'s horizon the stretch is a ten-thousandth of g. Coordinates that stay regular at the horizon (Gullstrand–Painlevé) picture space as a river flowing inward at the escape speed, which reaches c at the horizon (Hamilton & Lisle 2008). A spinning (Kerr) black hole adds an ergosphere, where nothing can stand still, and an inner horizon. **Per ground rule 1 this is a simulation**, and every input is measured: Sgr A* from H1, M87* from I5, Cygnus X-1 from I4, GW250114's remnant from I7, and the sky from real star catalogues and maps.

**Data.** No new observations. The background sky, all from the ESA Gaia archive: the 1.25 million Gaia DR3 stars with G < 11 as point stars (Hipparcos V and B−V replace Gaia's saturated photometry below V = 5, via `gaiadr3.hipparcos2_best_neighbour`), plus the G, BP and RP fluxes of all 1.8 billion fainter sources summed on the archive's server per 3.4′ HEALPix pixel (`GROUP BY source_id/2^35/4^2` over 192 sky chunks, ~10 min). A hips2fits image (H5's route) would count the point stars twice and comes with a non-linear stretch. Colours from the CIE 1931 colour-matching functions (CVRL) and Mamajek's colour–temperature table. Masses from H1, I4 and I5; I7 is not built yet, so GW250114's remnant comes from the LVK paper (PRL 135, 111403: 62.7 M☉, χ = 0.68 ± 0.01), Cyg X-1's spin from Zhao et al. 2021 (> 0.9985, contested).

**Method + visualizations.**
- Radial infall: proper time and Schwarzschild time against radius, from rest far away and from rest near the horizon; the faller's clock against a distant clock; the horizon-to-singularity time for each real black hole.
- Tidal acceleration across a 2 m body against radius for each black hole, with a ~10 g survival line → where spaghettification starts, outside or inside the horizon.
- What the faller sees: trace rays backward through Schwarzschild spacetime (I6's integrator) from the faller's local frame, with aberration and Doppler shift from their speed, to the real sky → a movie of the descent into Sgr A*: the shadow growing, the sky crowding into a bright, blueshifted ring ahead.
- What the friend sees: the faller's signals redshifting and fading exponentially.
- A Penrose (conformal) diagram and a river-model diagram with the faller's world line. Kerr extension: frame dragging in the ergosphere, and an orbit plunging from the innermost stable orbit.

**Check against.** Analytic results: proper time from horizon to singularity ≤ πGM/c³; the shadow edge at impact parameter √27 GM/c²; Hamilton & Lisle 2008's river model; Andrew Hamilton's published Black Hole Flight Simulator movies.

**Watch out for.** Aberration must be applied in the faller's local frame, not in Schwarzschild coordinates, or the view near the horizon is wrong. Rays skimming the photon sphere need adaptive steps. Real black holes are surrounded by hot gas (I5, I10); our clean view leaves it out, so say so. Inside the horizon Schwarzschild's coordinates swap roles; use Gullstrand–Painlevé or Kruskal coordinates for the interior. Found when this was built: the background is the Sun's sky, not the Galactic Centre's (no 3D catalogue gives the latter), so say so; a falling observer's dark region is only 42.1° in radius at the horizon, while a hovering one's fills the sky (aberration again); whatever the mass, a body falling from far away feels 10 g of stretch across 2 m exactly 0.095 s before the singularity, (2/3)√(L/a); the sky spans ~6 decades of brightness, so a log-like (asinh) tone curve is needed; the textbook arctan(U) Penrose compactification hides observers beyond a few r_s in the corners (arctan(|U|^½) does not); frame dragging is invisible in horizon-penetrating Kerr coordinates and shows in the Boyer–Lindquist angle against the distant stars.

**Builds on.** I6 (ray tracing), H1, I4, I5 and I7 (real masses), H5 (the sky), N6 (the conserved quantities that make the integration simple), I10 (tides on stars).

## N6. Killing vectors: the symmetries that make orbits solvable ★★☆

**Science question.** Why can we work out S2's orbit, a black hole's shadow and the innermost stable orbit with pencil and paper, when Einstein's equations are so hard? Because those spacetimes have symmetries, described by **Killing vector fields**, and every symmetry gives a conserved quantity. Can we find a metric's Killing vectors ourselves with computer algebra, use the conserved quantities to solve orbits, and test them on real data: S2, the EHT shadows, the cosmological redshift and black hole thermodynamics?

**Background.** A Killing vector field ξ points along a direction in which the metric does not change. It satisfies Killing's equation, `∇_μ ξ_ν + ∇_ν ξ_μ = 0`. Along any geodesic with four-velocity u, the quantity ξ·u stays constant: Noether's theorem in geometric form. Flat spacetime has 10 Killing vectors (translations, rotations and boosts). Schwarzschild spacetime has four: time translation ∂_t, which gives the conserved **energy per unit mass** E = (1 − 2GM/rc²) c² dt/dτ, and three rotations, which give the **angular momentum** L. With E and L an orbit reduces to motion in a one-dimensional effective potential. That gives the precession (N1), the photon sphere at r = 3GM/c², the shadow edge at impact parameter √27 GM/c² (I5, I6) and the innermost stable circular orbit at 6GM/c², where accreting gas has released 5.7% of its rest energy. A spinning Kerr black hole has only two Killing vectors, but also a hidden symmetry, a Killing *tensor*, whose conserved **Carter constant** (1968) keeps orbits regular rather than chaotic. In the ergosphere the time-translation Killing vector turns spacelike and E can be negative, which lets the Penrose process extract rotational energy, up to 29% of a maximally spinning black hole's mass. A horizon is a Killing horizon, and its surface gravity κ sets the Hawking temperature `T = ħκ / (2π c k_B)`, about 6 × 10⁻⁸ K × (M☉/M). An expanding FLRW universe has spatial Killing vectors but no timelike one, so photon energy is not conserved: that is the cosmological redshift, and it is why peculiar velocities decay as 1/a.

**Data.** Calculation first, with SymPy 1.14 (already in the shared venv): metric → Christoffel symbols → Killing equations. The real-data tests reuse earlier results: H1's S2 orbit and posterior samples, I5's ring diameter, I7's GW250114 remnant mass and spin, published X-ray-binary spins (I4), and D7's CMB temperature against redshift. *(Corrected 2026-10-19 when N6 was built: I7 did not exist yet, so N6 downloads LVK's own GW250114 posterior samples (NRSur7dq4, Zenodo 10.5281/zenodo.16877102, 27 MB); I4 measured masses, not spins, so the spins (19 systems, continuum fitting and reflection) were compiled from the papers and are strongly model-dependent.)*

**Method + visualizations.**
- Solve Killing's equation symbolically for flat spacetime (10 solutions), Schwarzschild (4), Kerr (2) and spatially flat FLRW (6, none timelike), and check each by substitution; also FLRW's conformal Killing vector a(t)∂_t.
- Effective potentials for particles and light in Schwarzschild and Kerr (prograde and retrograde); innermost stable orbit and photon orbit radii against spin, and the accretion efficiency from 5.7% to 42%.
- S2 as a test: along H1's fitted orbit, compute E and L at every observed epoch from the fitted positions and velocities, in Newtonian and GR-corrected form → constant within the errors; their values across the posterior samples. *(As built: a fitted orbit conserves E and L by construction, so N6 tests conservation season by season, fitting each season's positions (and spectra) alone, relative to the black hole as located by H1, with each team's frame drift: Kepler's second law and the energy, plus E, L and the pericentre shift over the posterior.)*
- Shadows from conserved quantities: the critical impact parameter for Sgr A* and M87* from H1's and I5's masses, against the EHT rings, with the small spread over spin.
- Thermodynamics: surface gravity, temperature, entropy (horizon area) and evaporation time for every black hole in the atlas, from GW250114's remnant to M87*; I7's area increase read as an entropy increase.
- Cosmology: why light redshifts in FLRW, and the T ∝ (1 + z) that D7 measured.

**Check against.** Standard results (Carroll 2004; Wald 1984): the Killing-vector counts, r_ISCO = 6GM/c², b = √27 GM/c², efficiencies of 5.7% and 42%, T_H = 6.2 × 10⁻⁸ K for 1 M☉; H1's measured orbit; I5's ring.

**Watch out for.** Real orbits are not exact geodesics: S2's velocities have errors, and H1's fitted reference-frame drift must be included, so "constant" means "consistent with constant". Kerr's Killing tensor does not come out of the vector equation; derive it separately or verify the known form. Symbolic simplification is slow for Kerr; substitute numbers early or use Boyer–Lindquist identities. Hawking radiation from astrophysical black holes is far colder than the CMB and has never been observed; it is a prediction, not a result. Counting is easier than solving: Taylor-expanding Killing's equation about one point in exact rational arithmetic gives a rigorous upper bound on the number of Killing vectors (and tensors), which meets the explicit solutions. A spacetime without Carter's tensor (Zipoy–Voorhees) still shows regular-looking orbits wherever orbits stay outside its barrier (KAM); the missing symmetry is proved by the count, not by the pictures.

**Builds on.** H1 (S2), I5 and I6 (shadows), I7 (the area theorem), I4 (spins), D7 (CMB temperature), N1 and N5.

## N7. Is Einstein still right? A scorecard of gravity tests ★☆☆

**Science question.** GR has been tested for over a century, from the laboratory to the whole universe, and has passed every time. How precisely, and where could it still fail? Can we compute GR's prediction for every classic test and plot each measurement against it?

**Background.** The tests fall into groups. The **equivalence principle**, that all bodies fall alike, has been tested from Eötvös (10⁻⁹) to the MICROSCOPE satellite: η(Ti, Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)] × 10⁻¹⁵ (Touboul et al. 2022). Its strong form says that gravitational binding energy falls like everything else. Lunar laser ranging tests it with the Earth and Moon falling toward the Sun, and the triple pulsar PSR J0337+1715, a neutron star and a white dwarf falling together toward an outer white dwarf, gives |Δ| < 2.1 × 10⁻⁶ (Voisin et al. 2020; 1.5–2.3 × 10⁻⁶ in Voisin et al. 2025, depending on the noise model). **Weak-field tests** measure γ and β (N1–N4). **Strong-field and dynamical tests** include binary pulsars (C9), gravitational waves (C1, I7; GW170817 showed in C4 that gravity travels at the speed of light to 1 part in 10¹⁵), S2's orbit (H1) and black hole shadows (I5). On the largest scales, dark matter and dark energy (sections D and K) mean either unseen ingredients or different gravity. Antimatter falls too: ALPHA-g measured antihydrogen's acceleration at 0.75 ± 0.13 ± 0.16 g (2023).

**Data.** Published results, each with its reference: Will's *Living Reviews in Relativity* review (2014) for the history, plus MICROSCOPE 2022, Voisin et al. 2020, Ciufolini et al. 2026, the Galileo redshift tests (2018), Cassini (2003), and this atlas's own measurements. APOLLO's public lunar laser ranging normal points (`tmurphy.physics.ucsd.edu/apollo/norm_pts.html`, 3,846 points 2006–2020) to illustrate what lunar ranging measures; modelling them yourself needs JPL DE440 and its lunar orientation kernel from NAIF (~45 MB), the IERS EOP 20 C04 series, and JPL's range-model memo (`ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf`), with the reflector coordinates in Park et al. 2021, Table 1. The triple pulsar's arrival times are public (Voisin et al. 2025, Zenodo 10.5281/zenodo.13899771). As built in 64-gravity-scorecard the download is ~0.2 GB with cross-checks (EPM2021 byte range from `ftp.iaaras.ru`, DE421/DE430, GPS velocities).

**Method + visualizations.**
- The equivalence-principle timeline: η from Eötvös (1908) through the torsion balances of Dicke, Braginsky and the Eöt-Wash group to MICROSCOPE, on a log scale.
- The γ timeline (with N2) and the β constraints, as in Will's review.
- The triple pulsar explained: why a violation would pull the inner orbit toward the outer white dwarf, and the size of that signal against Voisin et al.'s timing precision.
- Lunar ranging: round-trip times of ~2.5 s from APOLLO normal points, measured to millimetres, and the Nordtvedt signal GR forbids, a 13 m × η stretching of the Moon's orbit toward the Sun.
- One summary chart: every test in the atlas, its fractional precision against the strength of gravity it probes (GM/rc²), from the laboratory to black holes.

**Check against.** Will 2014's tables; the papers above.

**Watch out for.** Published bounds use different confidence levels (1σ, 68%, 95%) and conventions (η or Δ, γ − 1 or γ); convert before plotting. Many results are model-dependent; the PPN framework covers only a class of theories. Passing tests narrows the alternatives; it does not prove GR complete, and GR still cannot be reconciled with quantum mechanics. Will 2014's tables are not error-free: it quotes Lebach et al. (1995) at twice its published error, its α₁ limit from PSR J1738+0333 is the full width of a 95% interval, and its Fig. 5 draws pre-1950 probable errors as 1σ. Residual tests against a GR ephemeris (DE440, EPM2021) are consistency checks: the ephemerides were fitted to lunar ranges assuming GR. DE440's Moon sits 1–2 m ahead along its orbit of EPM2021, DE421 and DE430 during 2006–2020; a range model built from JPL's memo sees it as a 7-cm sin(anomaly) residual in APOLLO data, which EPM2021 does not show.

**Builds on.** N1–N4, C1, C4, C9, H1, I5, I7, M4 (the triple pulsar is a hierarchical triple).

---

# O. The expanding universe: Friedmann models and the Big Bang

In 1922 Alexander Friedmann solved Einstein's equations for a universe that is the same everywhere and in every direction, and found that it cannot stand still: it must expand or contract. Georges Lemaître (1927) connected that expansion with the redshifts of galaxies and later traced it back to a hot, dense beginning. The **Friedmann equation**

`H² = (ȧ/a)² = (8πG/3) ρ − k c²/a² + Λc²/3`

says how fast the scale factor a(t) grows given what the universe contains. Each ingredient dilutes differently as space expands: radiation as a⁻⁴, matter as a⁻³, curvature as a⁻², a cosmological constant not at all. In terms of today's density parameters, `H(a)² = H0² [Ω_r a⁻⁴ + Ω_m a⁻³ + Ω_k a⁻² + Ω_Λ]`.

Section D measured the expansion. This section works with the equations themselves and follows the universe from 380,000 years onward; section P goes back to the first second. With Planck 2018 parameters:
- ~51,000 years (z ≈ 3,400): matter overtakes radiation.
- ~372,000 years (z ≈ 1,090): atoms form and the CMB is released (O4).
- ~0.2–0.5 billion years (z ≈ 20–10): the first stars light up the 21-cm sky (O5).
- ~1 billion years (z ≈ 6): reionization ends (I15).
- ~7.7 billion years (z ≈ 0.6): the expansion starts to accelerate (D1); dark energy overtakes matter at ~10 billion years (z ≈ 0.3).
- 13.8 billion years: today.

**Suggested Friedmann-universe arc:** O1 (build a universe) → O2 (horizons) → O4 (recombination) → O5 (dark ages and cosmic dawn) → O6 (is the universe the same in every direction?) → O3 (why the night sky is dark).

## O1. Build a universe: solving the Friedmann equations ★☆☆

**Science question.** Give a universe some matter, radiation, curvature and a cosmological constant: how does it evolve? Some recollapse into a Big Crunch, some expand forever, and some never had a Big Bang at all. Can we solve the Friedmann equations for all of them, map which universes our measurements allow, and find our own universe's age, past and future?

**Background.** From the Friedmann equation and the fluid equation, `ρ̇ = −3H(ρ + p/c²)`, each component with pressure p = wρc² dilutes as a^(−3(1+w)). The acceleration equation, `ä/a = −(4πG/3)(ρ + 3p/c²) + Λc²/3`, shows that ordinary matter slows the expansion and Λ speeds it up. Exact solutions: a ∝ t^(1/2) for radiation, a ∝ t^(2/3) for matter (Einstein–de Sitter), a ∝ t for an empty universe (Milne), a ∝ exp(Ht) for Λ alone (de Sitter); a closed matter universe follows a cycloid into a Big Crunch. Einstein's 1917 static universe balanced matter against Λ exactly, and Eddington showed in 1930 that it is unstable. Lemaître's "loitering" models hover near that balance. With Λ large enough there is no Big Bang, only a bounce, which the existence of quasars at high redshift rules out. D1 measured Ω_m = 0.295 ± 0.054 and Ω_Λ = 0.609 ± 0.081 from supernovae alone; adding the CMB and BAO gives Ω_k = 0.0007 ± 0.0019 (Planck 2018).

**Data.** D1's Pantheon+ likelihood with curvature free, D5's DESI DR2 BAO table and compressed CMB prior, D3's CAMB fit, and B6's oldest-star and B10's white dwarf ages as lower limits on the age. DESI's released ΛCDM + Ω_k chains and best fits (data.desi.lbl.gov `bao-cosmo-params`, `base_omegak`, ~0.14 GB) serve as the answer key; the rest (~60 MB) is reused from D1, D5, D8 and I3.

**Method + visualizations.**
- Integrate a(t) numerically over a grid of (Ω_m, Ω_Λ, Ω_r): a gallery of open, closed, flat, loitering, bouncing and recollapsing universes, all scaled to today's H0.
- The classic Ω_m–Ω_Λ plane: the flat line, the boundaries between acceleration and deceleration and between recollapse and eternal expansion, and the "no Big Bang" corner, with D1's supernova contours, the BAO and CMB contours, and where they cross.
- The age of the universe across the plane, with the oldest stars as a lower bound: the "age crisis" of the 1990s, when a matter-only universe was younger than its stars.
- An interactive "build-a-universe" page: sliders for Ω_m, Ω_Λ, Ω_r and w, showing a(t), H(z), the age and the fate, with our universe marked.
- The equality epochs: when radiation gave way to matter, and matter to dark energy.

**Check against.** Planck 2018: age 13.787 ± 0.020 Gyr with BAO (13.797 ± 0.023 for TT,TE,EE+lowE+lensing alone), Ω_k = 0.0007 ± 0.0019 (+lensing+BAO); DESI DR2 Table V: Ω_k = 0.025 ± 0.041 (BAO alone), 0.0023 ± 0.0011 (with the CMB); D1's contours; analytic ages (2/(3H0) for Einstein–de Sitter, 1/H0 for Milne).

**Watch out for.** The "no Big Bang" and loitering boundaries need radiation included near a → 0. Ages scale with 1/H0, so the Hubble tension (D1 against D3) moves them by up to a billion years. Keep the interactive page small and working offline.

**Builds on.** D1 (Ω_m–Ω_Λ), D3 (CAMB), D5 (BAO), D8 (H(z) and the future), B6 and B10 (ages).

## O2. The cosmic light cone: horizons and faster-than-light recession ★☆☆

**Science question.** How can the observable universe be 46 billion light-years in radius if it is only 13.8 billion years old? Why can we see galaxies that are receding faster than light? And which galaxies could we never reach, even travelling at light speed? Can we redraw the spacetime diagrams that answer these questions (Davis & Lineweaver 2004) with today's cosmology, and put real objects on them?

**Background.** In an expanding universe "distance" has several meanings: the proper distance now, the proper distance when the light left, the luminosity distance (D1) and the angular-diameter distance. The **Hubble sphere**, where the recession speed equals c, is now 14.5 billion light-years away (c/H0). Galaxies beyond it (z > 1.48) recede faster than light, which special relativity does not forbid because space itself is expanding. We still see them, because the Hubble sphere grows, so light that was initially carried away from us eventually approaches. The **particle horizon**, the farthest that light could have travelled since the Big Bang, is now 46 billion light-years away. The **event horizon**, beyond which light we emit today can never arrive, is 16.6 billion light-years away: galaxies now beyond z ≈ 1.85 are forever out of reach. Angular sizes have a twist: an object of fixed size looks smallest at z ≈ 1.6 and larger beyond, because its light left when the universe was small.

**Data.** D3's Planck parameters (or `astropy.cosmology.Planck18`). Real objects from earlier projects: D2's SDSS galaxies, I3's quasars, D1's supernovae, C5 and I9's gravitational-wave events, J8's JWST galaxies at z ≈ 14, and the CMB. D5's BAO angular-diameter distances, to show the turnover in data. (As built: C5, I9 and J8 did not exist yet, so the GWOSC event list and the JADES DR4 redshift table are read directly; the latter is the first FITS extension of a 94 MB file, 4.4 MB by HTTP byte range. Davis & Lineweaver's arXiv source carries their figures as vector EPS, which can be read back for an exact check. Download ~5.5 MB when D1, D2, D5 and I3 have their caches, ~85 MB otherwise.)

**Method + visualizations.**
- Davis & Lineweaver's three-panel spacetime diagram, recomputed: proper distance, comoving distance and conformal time against cosmic time, with our past light cone (teardrop-shaped in proper distance), the Hubble sphere, the particle and event horizons, and the world lines of galaxies at z = 1, 3, 10 and 1,090.
- Real objects placed where their light was emitted: the CMB, J8's galaxies, I3's quasars, supernovae and gravitational-wave events.
- Each object's recession speed now and when its light left, against c.
- The angular size of a fixed 1 Mpc ruler against redshift, with DESI's BAO measurements on it.
- Reachable or not: the share of the observable universe we could still send a signal to (about 5% of its volume).

**Check against.** Davis & Lineweaver 2004, Figure 1. With Planck 2018: particle horizon 46.2 billion light-years, event horizon 16.6, recession speed c at z = 1.48, smallest angular size at z = 1.59, event horizon at z = 1.85.

**Watch out for.** Mixing up distance measures is the whole topic, so label every axis with which distance it is. "Faster than light" means proper distance growing faster than c in cosmic time; nothing overtakes a light beam locally. The numbers depend on the cosmology, so quote the parameters used. Davis & Lineweaver's Figure 1 has no radiation (with it, their particle horizon would be 45.3 rather than 46.2 billion light-years). The special-relativistic and v = cz curves of their Figure 5 use D_L = D, not the D_L = D(1+z) their text describes, and their 23σ follows from the drawn curves.

**Builds on.** D1, D3 and D5 (distances), D7 (the expansion is real), J8 (the most distant galaxies), I3 (quasars).

## O3. Why is the night sky dark? Olbers' paradox and the cosmic optical background ★★☆

**Science question.** If the universe were infinite, eternal and unchanging, every line of sight would end on a star, and the night sky would blaze as bright as the Sun's surface. It is dark. Why? And how dark is it really, measured from beyond Pluto, where the glow of the inner Solar System's dust fades away?

**Background.** Olbers' paradox (1823, already discussed by Kepler) is resolved mainly by the universe's finite age: light from beyond ~14 billion light-years has not had time to arrive, and stars shine for only a limited time, so there is nowhere near enough starlight to fill the sky (Harrison, *Darkness at Night*, 1987). Expansion adds redshift and dilution. All the light ever emitted by stars and accreting black holes, redshifted and accumulated, makes the **extragalactic background light**, a few tens of nW m⁻² sr⁻¹ split between the optical and the far infrared (light absorbed and re-emitted by dust). New Horizons' LORRI camera, far from the zodiacal light, measured the cosmic optical background at 11.08 ± 1.65 nW m⁻² sr⁻¹, consistent with the integrated light of known galaxies, 8.17 ± 1.18 (Postman et al. 2024), after earlier hints of an excess (Lauer et al. 2021, 2022). Gamma rays from distant blazars are absorbed by this light, which gives an independent measurement (Fermi-LAT Collaboration 2018). For comparison, the CMB's intensity is ~1,000 nW m⁻² sr⁻¹: most of the universe's light is the Big Bang's afterglow.

**Data.** Galaxy number counts from the far ultraviolet to the far infrared, with the integrated galaxy light at each wavelength (Driver et al. 2016, published tables). The New Horizons measurements (Lauer et al. 2021, 2022; Postman et al. 2024). The COBE DIRBE/FIRAS infrared background (Hauser & Dwek 2001 review) and the Fermi-LAT gamma-ray constraints (2018). Optionally, our own counts from J9's Hubble Ultra Deep Field catalogue.

**Method + visualizations.**
- The paradox in numbers: sky brightness for an infinite static universe (the surface brightness of a star), for one of finite age (stellar emission integrated over a light-travel shell), and for an expanding one (with D7's (1 + z)⁻⁴ surface-brightness dimming).
- Integrated galaxy light: counts at each wavelength integrated over magnitude → the total, which converges at faint magnitudes, so most of the light is in galaxies already seen.
- The extragalactic background spectrum from ultraviolet to far infrared, with the New Horizons point, the galaxy light and the gamma-ray constraints, and the CMB beside it for scale.
- How dark is dark: the sky from Earth (airglow, zodiacal light), from beyond Pluto, and the cosmic floor.

**Check against.** Postman et al. 2024: 11.08 ± 1.65 (background) and 8.17 ± 1.18 (galaxy light) nW m⁻² sr⁻¹; Driver et al. 2016's integrated galaxy light.

**Watch out for.** Counts must be corrected for incompleteness at faint magnitudes and for the extended light of galaxies, and deep fields differ by ~10–20% (cosmic variance). Absolute sky measurements fight scattered starlight and dust light from our own Galaxy, and the New Horizons values depend on those models. Convert units carefully between AB magnitudes, janskys and nW m⁻² sr⁻¹.

**Builds on.** D7 (surface-brightness dimming), J9 (the deep field), J3 (cosmic star formation), D3 (the CMB).

## O4. Recombination: the moment the universe became transparent ★★☆

**Science question.** For its first 380,000 years the universe was a glowing fog of plasma, as opaque as the inside of the Sun. Then electrons and protons combined into atoms, and light streamed free: the CMB. Why did this happen at ~3,000 K, when hydrogen's ionization energy corresponds to 158,000 K? Can we compute the ionization history ourselves, from the simplest equation to the full physics, and match the redshift of last scattering that Planck measured?

**Background.** In equilibrium, the Saha equation gives the ionized fraction. There are ~1.6 billion photons per baryon, so even at a few thousand kelvin the high-energy tail of the blackbody holds enough photons to keep hydrogen ionized; by Saha, recombination is half complete only at ~3,760 K (z ≈ 1,380). Real recombination is slower still. A direct recombination to the ground state emits a photon that ionizes another atom, and Lyman-α photons are reabsorbed at once, so atoms reach the ground state only through the slow two-photon decay of the 2s level (8.2 per second) or by Lyman-α photons redshifting out of the line (Peebles 1968; Zel'dovich, Kurt & Sunyaev 1968). The **visibility function**, the probability that a CMB photon last scattered at redshift z, peaks at z* ≈ 1,090. Helium recombined earlier, in two steps (z ≈ 6,000 and 2,500). Much later the first stars and galaxies reionized the gas; Planck's Thomson optical depth τ = 0.054 ± 0.007 puts the midpoint at z ≈ 7.7 (I15).

**Data.** Mostly calculation. CAMB (D3, in the shared venv) for the full recombination history (RECFAST). Planck 2018 parameters: z* = 1089.92 ± 0.25, z_drag = 1059.94 ± 0.30, τ = 0.054 ± 0.007 (paper tables, or the chains in the Planck Legacy Archive). FIRAS's T₀ = 2.7255 K and D3's own Ω_b h².

**Method + visualizations.**
- Saha ionization fractions against redshift for hydrogen and helium.
- Peebles' three-level atom as one differential equation, solved with SciPy, against Saha and against CAMB → recombination "freezes out", leaving a residual ionization of a few × 10⁻⁴.
- The optical depth and visibility function → z*, and how it shifts with Ω_b h² and T₀.
- The age at last scattering, and the last-scattering shell on O2's light cone.
- The full history x_e(z) including reionization (CAMB's tanh model) for Planck's τ, with I15's quasar constraints.

**Check against.** Planck 2018: z* = 1089.92 ± 0.25, z_drag = 1059.94 ± 0.30; age at z* ≈ 372,000 years; Saha's halfway point at ~3,700 K (z ≈ 1,370; Ryden, *Introduction to Cosmology*).

**Watch out for.** The three-level atom is good to a few percent near z*; Planck needs many-level atoms and radiative transfer (HyRec, CosmoRec). Helium must be included to get x_e right above z ≈ 1,500. Don't confuse z* (when photons last scattered) with z_drag (when baryons were released, which sets D5's BAO ruler).

**Builds on.** D3 (CMB, CAMB), K1 (the baryon density), D5 (the sound horizon), I15 (reionization), O2 (the light cone).

## O5. The dark ages and cosmic dawn: hydrogen's 21-cm signal ★★☆

**Science question.** Between the CMB's release and the first stars there was nothing to see, only cooling hydrogen: the cosmic dark ages. The first stars then flooded the gas with ultraviolet light, which made it absorb CMB radiation at hydrogen's 21-cm line, redshifted today to ~50–100 MHz. In 2018 the EDGES experiment reported that dip at 78 MHz, more than twice as deep as standard physics allows. In 2022 SARAS 3 looked and did not find it. Can we compute what standard cosmology predicts, fit the EDGES data ourselves, and see how fragile such a measurement is?

**Background.** The 21-cm line comes from a flip of the hydrogen electron's spin relative to the proton's. Its excitation is described by a spin temperature T_s: gas with T_s below the CMB temperature absorbs CMB photons, gas above it emits. Against the CMB the signal is roughly `δT_b ≈ 27 mK × x_HI × (1 − T_CMB/T_s) × [(1 + z)/10]^(1/2)` for Planck parameters. After recombination, leftover free electrons kept the gas at the CMB's temperature until z ≈ 150; then the gas cooled faster, as (1 + z)², than the radiation, as (1 + z). Ultraviolet light from the first stars tied T_s to the cold gas (the Wouthuysen–Field effect) at z ≈ 15–20, producing an absorption dip of at most ~0.2 K in standard models, before X-rays heated the gas. EDGES found 0.5 (+0.5/−0.2) K centred at 78 ± 1 MHz (z ≈ 17; Bowman et al. 2018), which would need colder gas (perhaps cooled by dark matter) or an extra radio background. Others argued the dip is an artefact of the fitting (Hills et al. 2018), and SARAS 3 rejected the EDGES profile at 95.3% confidence (Singh et al. 2022). Interferometers such as HERA instead measure the signal's fluctuations, and have shown the gas was heated by z ≈ 8 (HERA Collaboration 2023).

**Data.** EDGES data releases (`loco.lab.asu.edu/edges/edges-data-release/`): Data Release 1 has the reduced, averaged spectra and fitted profile of Bowman et al. 2018; Data Release 3 the raw data as reanalysed by Murray et al. (2026). Published SARAS 3 and HERA limits. D3's Planck parameters and O4's residual ionization for the thermal history.

**Method + visualizations.**
- The thermal history: CMB, gas and spin temperatures against redshift (and against observed frequency, 1420/(1 + z) MHz), from Compton heating by O4's residual electrons plus a simple model of the first stars' ultraviolet and X-ray light.
- The predicted global signal δT_b(ν), with a band over plausible first-star models and the standard-physics floor.
- Refit EDGES Data Release 1: their five-term foreground model plus a flattened Gaussian → the 78 MHz dip; then other foreground models and a sinusoidal systematic, and how the profile changes (Hills et al. 2018).
- EDGES, SARAS 3 and the standard prediction on one plot, with HERA's lower limit on the gas temperature.

**Check against.** Bowman et al. 2018: 78 ± 1 MHz, 0.5 (+0.5/−0.2) K, width 19 (+4/−2) MHz; Singh et al. 2022: rejected at 95.3%.

**Watch out for.** The sky's synchrotron emission is ~10⁴ times brighter than the signal, and any frequency dependence in the antenna's beam turns smooth foregrounds into fake features. The ionosphere and ground plane matter too. A fit with enough free foreground terms can hide a dip or create one, so report how the result depends on the model. The thermal history needs O4's residual ionization, not zero.

**Builds on.** O4 (residual ionization), H2 (the 21-cm line), D3 (cosmology), J8 (the first galaxies), I15 (reionization).

## O6. Is the universe the same in every direction? The cosmic dipole ★★☆

**Science question.** Modern cosmology assumes that, averaged over large scales, the universe looks the same in every direction: the cosmological principle. The CMB is hotter in one direction by 1 part in 1,000 (3.36 mK), which we read as our own motion at 370 km/s. If that is right, distant quasars and radio galaxies should show the same motion, as a slight excess of sources in that direction. They do show an excess, but about twice as large as expected. Is the cosmological principle in trouble?

**Background.** Moving through a sea of sources, aberration crowds them toward the direction of motion and Doppler boosting brightens them there, so a survey with a flux limit counts more sources ahead. Ellis & Baldwin (1984) worked out the expected dipole amplitude, `D = [2 + x(1 + α)] v/c`, where x is the slope of the integral source counts and α the spectral index; for typical quasars and v = 370 km/s it is under 1%. Planck measured the CMB dipole as 369.82 ± 0.11 km/s toward Galactic coordinates (l, b) = (264.0°, 48.3°). Secrest et al. (2021) counted over a million mid-infrared quasars from CatWISE2020 and found a dipole within ~27° of the CMB's direction but twice as large, rejecting the expected amplitude at 4.9σ. Radio samples point the same way, and a 2025 review quotes over 5σ in several independent samples (Secrest et al. 2025, arXiv:2505.23526). A 2025 reassessment that models clustering and the survey footprint in simulations finds 3.6σ for CatWISE (Bashir et al. 2025). If the excess is real, matter and the CMB do not share one rest frame.

**Data.** The Secrest et al. (2022) release on Zenodo (record 6784602, ~2.4 GB of archives with the CatWISE2020 quasar and NVSS radio samples and their masks). The NVSS catalogue (VizieR `VIII/65`). Planck's dipole (D3). Optionally Gaia's quasar candidates (H4) as a third sample.

**Method + visualizations.**
- Sky maps of quasar counts (HEALPix, `healpy`), with the Galactic plane masked and the CMB dipole direction marked.
- Measure the dipole by least squares and by Poisson maximum likelihood, correcting for the way WISE's depth varies with ecliptic latitude → amplitude and direction.
- The expected kinematic dipole from the catalogue's own x and α for v = 370 km/s → the ratio of measured to expected.
- Significance by simulation: mock skies with the kinematic dipole, shot noise and clustering (lognormal), seen through the same mask → how often a dipole this large appears by chance.
- Repeat for NVSS, and for Gaia's quasars.

**Check against.** Secrest et al. 2021: twice the kinematic expectation, 4.9σ; Bashir et al. 2025: 3.6σ with clustering included; Planck: 369.82 ± 0.11 km/s.

**Watch out for.** Anything that varies smoothly across the sky can mimic a dipole: dust extinction, WISE's scanning depth, star contamination, the Galactic mask. Clustering of nearby sources adds a real but non-kinematic dipole. x and α must be measured at the flux limit. This is an open question; present it as one.

**Builds on.** D3 (the CMB dipole), H4 (Gaia's quasars), I3 (quasars), I18 (radio surveys), J4 (our motion relative to nearby galaxies).

---

# P. The first second: inflation and particles after the Big Bang

Run the expansion backward and the universe gets hotter, with temperature rising as 1/a. At high enough temperatures atoms break apart (O4), then nuclei (K1), then protons and neutrons dissolve into quarks and gluons. Radiation dominates, and the Friedmann equation becomes a simple clock: `t ≈ 2.4 s × g*^(−1/2) × (1 MeV / kT)²`, where g* counts the particle species that are relativistic at temperature T (1 MeV corresponds to 1.16 × 10¹⁰ K). Much of particle physics played out in the first second. This section follows it with laboratory and sky data.

- **Inflation** (time unknown, perhaps ~10⁻³⁶–10⁻³² s): a burst of accelerated expansion makes the universe flat and uniform and seeds its structure (P1). **Reheating** turns its energy into hot particles.
- **Baryogenesis** (epoch unknown): a tiny excess of matter over antimatter appears, one extra quark per ~270 million quark–antiquark pairs (P5).
- **Electroweak crossover**, ~10⁻¹¹ s, ~160 GeV: the Higgs field takes its present value, and W and Z bosons, quarks and leptons acquire mass (P7).
- **Dark matter freeze-out**, if dark matter is a WIMP, at T ≈ m/20 (P4).
- **Quark–hadron crossover**, ~10⁻⁵ s, ~155 MeV: quarks and gluons bind into protons, neutrons and pions (P6).
- **Neutrino decoupling**, ~1 s, ~1 MeV: neutrinos stop interacting and stream freely to this day (P3).
- **Electron–positron annihilation**, a few seconds: heats the photons relative to the neutrinos.
- **Nucleosynthesis**, ~3–20 minutes: deuterium and helium form (K1). Then section O takes over.

P2 turns this timeline into numbers.

**Suggested early-universe arc:** P1 (inflation) → P2 (the thermal history) → P7 (the Higgs) → P6 (quarks to hadrons) → P5 (antimatter) → P3 (relic neutrinos) → P4 (dark matter freeze-out) → K1 (nuclei) → O4 (atoms) → O5 (the first stars).

## P1. Inflation on trial ★★☆

**Science question.** Why is the universe so flat and so uniform, with CMB regions that were never in contact matching in temperature to 1 part in 100,000? Inflation's answer is that a brief burst of exponential expansion stretched a tiny patch to cosmic size, and that its quantum jitters became the seeds of galaxies. Different models of the field that drove it leave different patterns in the CMB. Which models survive Planck, BICEP/Keck, ACT and SPT, and which did the 2025 data put under pressure?

**Background.** Without inflation, the CMB sky splits into thousands of patches a degree or two across that could never have exchanged light (the horizon problem), and the density must have matched the critical value to ~1 part in 10⁶⁰ at the Planck time (the flatness problem). In slow-roll inflation a scalar field φ rolls slowly down a potential V(φ) while space expands by N ≈ 50–60 e-folds. Quantum fluctuations become density perturbations with a nearly scale-invariant spectrum, `P(k) ∝ k^(n_s − 1)`, and gravitational waves whose relative strength is the tensor-to-scalar ratio r. In terms of the slow-roll parameters ε = (M_Pl²/2)(V′/V)² and η = M_Pl² V″/V, n_s ≈ 1 − 6ε + 2η and r ≈ 16ε. At N = 60: m²φ² gives n_s = 0.967, r = 0.13; λφ⁴ gives 0.95 and 0.26 (16/N = 0.27 is the cruder estimate); Starobinsky's R² model (and Higgs inflation) gives 0.967 and 0.003. Planck 2018 measured n_s = 0.9649 ± 0.0042, and BICEP/Keck (BK18) limits r < 0.036 (95%), which excludes φ² and φ⁴. ACT DR6 combined with Planck, CMB lensing and DESI DR1 BAO shifted n_s up to 0.9743 ± 0.0034 (Louis et al. 2025), putting Starobinsky inflation near the 2σ edge; SPT-3G's CMB-SPA combination (SPT + Planck + ACT) finds 0.9679 ± 0.0033, and 0.9726 ± 0.0028 with DESI DR2 (values after the 2026 bug fix; SPT-3G alone gives 0.951 ± 0.011). Primordial gravitational waves would leave a swirling "B-mode" pattern in the CMB's polarization, the target of BICEP/Keck, the Simons Observatory and LiteBIRD.

**Data.** MCMC chains: BICEP/Keck BK18 data products (`bicepkeck.org`, October 2021 release, with bandpowers, likelihood and chains; LAMBDA's mirror `lambda.gsfc.nasa.gov/data/suborbital/BICEPK_2021/` is faster, and BK18's own chains are 66.7 GB, or 2.7 GB without data files, but not needed because ACT's tensor runs all include BK18), ACT DR6 chains on NASA's LAMBDA (`lambda.gsfc.nasa.gov/product/act/act_dr6.02/`), and Planck 2018 chains from the Planck Legacy Archive (the 11.8 GB grid zip on IRSA's mirror; single runs can be read out of it by HTTP range requests). SPT-3G D1 chains on LAMBDA (ΛCDM only, no r; read the same way out of 6.2 GB and 0.8 GB zips, the CMB-SPA runs from the 2026 bug-fix release). BICEP/Keck's 2026 component-separated maps (BK-XX) are not public, only figures; BK18's released per-bin CMB/dust/synchrotron decomposition shows where B-modes hide instead. D3's own n_s. Download ~0.23 GB.

**Method + visualizations.**
- The problems inflation solves: the comoving Hubble radius shrinking during inflation and growing afterwards, and the number of causally disconnected CMB patches without it.
- Solve the inflaton's equation of motion with the Friedmann equation for each potential (φ², φ⁴, φ^(2/3), natural, Starobinsky, hilltop, α-attractors) → the end of inflation and N; exact n_s and r from the Mukhanov–Sasaki equation for the perturbations, against the slow-roll formulas. (A labelled calculation.)
- The keeper: the n_s–r plane with Planck, BK18, ACT DR6 and SPT-3G contours from the chains, and each model as a band for N = 50–60.
- What r means: the energy scale of inflation, V^(1/4) ≈ 1.02 × 10¹⁶ GeV × (r/0.01)^(1/4), and the field excursion (the Lyth bound).
- B-modes: the predicted BB spectrum for r = 0.036 and 0.003 against lensing B-modes and the dust foreground, with BK18's bandpowers.
- Our own n_s: where D3's TT-only fit falls on the plane.

**Check against.** Planck 2018 (inflation paper): n_s = 0.9649 ± 0.0042; BK18: r < 0.036 (95%), and r < 0.035 with Planck + BAO; ACT DR6 (P-ACT-LB): n_s = 0.9743 ± 0.0034, and r < 0.038 with BK18 (Calabrese et al. 2025); the analytic slow-roll predictions.

**Watch out for.** n_s has shifted between data combinations by about 2σ, partly through how each handles the largest scales and the optical depth τ, so do not treat one combination as final. N depends on the unknown reheating history, so each model is a band, not a point. The B-mode limit hinges on modelling the dust foreground. Inflation is the leading explanation, not a proven one; alternatives such as bouncing cosmologies are still discussed. BK18's bandpower window functions include the maps' filtering suppression, so the likelihood's bandpowers sit ~7× below the released display bandpowers. ACT's and SPT's Cobaya chains include burn-in. LAMBDA's ACT chains use the "ACT-lite" likelihood (P-ACT-LB n_s 0.0005 below the paper's). SPT-3G never fitted r, so its n_s–r contours need an approximation (n_s posterior × BK18's r likelihood: 6% on r₉₅, checked on ACT).

**Builds on.** D3 (n_s, CAMB), O1 and O2 (flatness and horizons), P2 (reheating into a hot universe), I8 (gravitational-wave backgrounds), D5 (DESI).

## P2. From quarks to atoms: the thermal history of the universe ★☆☆

**Science question.** What was the universe made of at each moment of its first minutes, and how fast was it cooling? Can we count the Standard Model's particles species by species, add what lattice QCD says about the quark–hadron transition, and turn the Friedmann equation into a calendar of the first second?

**Background.** In the radiation era the energy density is `ρ = (π²/30) g* T⁴` (in natural units), where g* counts relativistic degrees of freedom: 2 for the photon, 7/8 for each fermion state. Above the electroweak scale the full Standard Model gives g* = 106.75. As the universe cools below a particle's mass, that particle annihilates and hands its entropy to what remains, so g* falls in steps: top quarks, the Higgs, W and Z, bottom and charm quarks, tau leptons. Then comes the big drop at the quark–hadron crossover, from ~62 (quarks, gluons and leptons) to ~17 (pions and leptons), at T_pc = 156.5 ± 1.5 MeV in lattice QCD (HotQCD 2019). After e± annihilation only photons and the now colder neutrinos remain (g* = 3.36). The Friedmann equation gives `t ≈ 2.42 s × g*^(−1/2) × (T / 1 MeV)^(−2)`. Neutrinos decouple at ~1 MeV, just before e± annihilation heats the photons, so the neutrino background is colder: T_ν/T_γ = (4/11)^(1/3) (P3). In the Standard Model the electroweak transition is a smooth crossover at 159.5 ± 1.5 GeV (D'Onofrio & Rummukainen 2016), not a sharp phase transition, which matters for P5.

**Data.** The Saikawa & Shirai (2018) tables of g*ρ(T), g*s(T) and the equation of state w(T), including lattice QCD (`member.ipmu.jp/satoshi.shirai/EOS2018`). Particle masses and spins from the PDG via the `particle` package (pure Python). Lattice crossover temperatures (HotQCD 2019; Borsányi et al. 2020) and the electroweak crossover (D'Onofrio & Rummukainen 2016). FIRAS's T₀ and K1's nucleosynthesis epoch for the end of the story. Also used (all small; checked 2026-10-20): Laine, Meyer & Schröder's 2015 Standard Model tables (`www.laine.itp.unibe.ch/eos15/`, http only), Borsányi et al.'s 2016 g*(T) spline table and 2+1+1 trace-anomaly fit (arXiv:1606.07494 source), the HotQCD 2014 (Table I) and Borsányi et al. 2014 lattice tables (arXiv sources), Escudero's 2020 decoupling table (arXiv:2001.04466 ancillary file) and the electron-mass rate coefficients of Escudero, Jackson, Laine & Sandner 2025 (nudec_BSM v2 on GitHub). About 3 MB in all.

**Method + visualizations.**
- Build g*(T) ourselves from the PDG particle list (ideal Fermi and Bose gases with their masses, integrated numerically) from 10 TeV down to 10 keV; compare with Saikawa & Shirai, who include interactions and lattice QCD → the difference is the QCD crossover.
- The staircase plot: g*ρ and g*s against T, each step labelled with the particle that disappears.
- The cosmic calendar: t(T) from the electroweak crossover (~10⁻¹¹ s) through the quark–hadron crossover (~10⁻⁵ s), neutrino decoupling (~1 s), e± annihilation and nucleosynthesis (minutes), as a log–log poster with the energies of the LHC, RHIC and the Sun's core marked.
- Entropy conservation through e± annihilation → T_ν/T_γ = (4/11)^(1/3) and today's 1.945 K neutrino background.
- The horizon size and the mass inside it at each epoch (a few solar masses at the QCD epoch: 3.7 M☉ at 156.5 MeV and 2.6 M☉ at the softest point, w = 0.229 at 178 MeV; Byrnes et al. 2018 put the resulting primordial-black-hole peak near 0.7 M☉; relevant for K4).

**Check against.** g* = 106.75 at high temperature, 10.75 at a few MeV, and 3.36 (energy) and 3.91 (entropy) today; HotQCD 2019: T_pc = 156.5 ± 1.5 MeV. These are ideal-gas values: with interactions Saikawa & Shirai give 104.3 ± 1.0 at 10 TeV, and with non-instantaneous neutrino decoupling 3.383 and 3.931 today (N_eff = 3.044; Escudero 2020).

**Watch out for.** g*ρ and g*s differ once species have different temperatures (after neutrino decoupling); use the right one in each formula. The ideal-gas count fails near the QCD crossover, where the strongly interacting plasma needs lattice QCD. Neutrino decoupling is not instantaneous (N_eff = 3.044, not 3; P3). Perturbative QCD converges poorly below ~10 GeV: its unknown g⁶ constant has to be matched to lattice data (Borsányi et al. 2016: q_c = −3000). Saikawa & Shirai's tabulated curve runs up to 0.5% below Escudero's decoupling run between 3 MeV and 0.2 MeV (beyond its quoted 0.1%) but agrees today. Hadron-gas lists differ on the broad f0(500) and K0*(700) (~1% at 100–150 MeV).

**Builds on.** K1 (nucleosynthesis), O1 (the Friedmann equation), O4 (recombination), K4 (primordial black holes), and P3, P6, P7.

## P3. The cosmic neutrino background ★★☆

**Science question.** One second after the Big Bang, neutrinos stopped interacting, and they have streamed through the universe ever since: 336 per cubic centimetre at 1.95 K. No one has detected them directly, but the CMB and the clustering of galaxies feel them. Can we see that, count the neutrino species from the CMB, and understand why cosmology's upper limit on neutrino mass now presses against the minimum that laboratory oscillations require?

**Background.** Relic neutrinos add to the radiation density in the early universe, measured as the effective number of neutrino species N_eff. The Standard Model predicts N_eff = 3.044, slightly above 3 because neutrinos had not fully decoupled when electrons and positrons annihilated (Bennett et al. 2021). More radiation speeds up the early expansion, which changes the CMB's damping tail and BBN's helium (K1). Planck finds N_eff = 2.99 ± 0.17; ACT DR6 with Planck, lensing and DESI, 2.86 ± 0.13. Neutrinos have mass: oscillation experiments measure differences of squared masses (Δm²₂₁ ≈ 7.4 × 10⁻⁵ eV², |Δm²₃₁| ≈ 2.5 × 10⁻³ eV²; NuFIT 6.0), which set a minimum total mass of ~0.059 eV (normal ordering) or ~0.10 eV (inverted). Massive neutrinos stream out of small structures and suppress clustering. DESI DR2 BAO plus the CMB limit Σm_ν < 0.064 eV (95%, in ΛCDM), excluding the inverted ordering's minimum and approaching the normal ordering's. KATRIN's direct laboratory limit is 0.45 eV (90%, 2025). A cosmological limit below the oscillation minimum would point to new physics or to a problem with the cosmological model.

**Data.** NuFIT 6.0's one-dimensional Δχ² tables (`www.nu-fit.org`, http only). Planck 2018 likelihood data as released for Cobaya's native likelihoods (GitHub `CobayaSampler/planck_native_data`: plik_lite, low-ℓ TT and EE; `planck_supp_data_and_covmats`: lensing), ACT DR6's lensing likelihood (LAMBDA; a 361 MB archive whose 2.9 GB of correction matrices reduce to 2 MB), BAO 2018 and DESI DR2 BAO (`CobayaSampler/bao_data`), DES-SN5YR (D5's copy). Answer keys: DESI DR2's released chains and best fits (`data.desi.lbl.gov`), Planck 2018's parameter-grid chains (ranged reads from IRSA's 11.8 GB zip), ACT's P-ACT-LB chain (LAMBDA). D3's CAMB setup, KATRIN's published limit, P2's thermal history. ~0.95 GB in all.

**Method + visualizations.**
- The relic neutrinos today: number density and temperature from P2's entropy argument, and their typical speeds today for masses of 0.01–0.1 eV.
- N_eff in the CMB: CAMB spectra for N_eff = 2, 3.044 and 4 at fixed acoustic scale, showing the damping-tail shift; then D3's TT spectrum refitted with N_eff free → our own constraint, against Planck's and ACT's.
- The neutrino mass plot: Σm_ν against the lightest neutrino mass for both orderings from NuFIT, with the cosmological limits (Planck, Planck + DESI DR1, DR2) and KATRIN.
- Our own bound: D3's CMB fit plus D5's DESI DR2 BAO with Σm_ν free (CAMB) → how much BAO tightens it, and how much evolving dark energy (D5's w0wa) loosens it.
- Matter power spectra with and without massive neutrinos: the small-scale suppression, about 8 Ω_ν/Ω_m.

**Check against.** Bennett et al. 2021: N_eff = 3.0440; Planck 2018: 2.99 ± 0.17; ACT DR6: 2.86 ± 0.13; DESI DR2: Σm_ν < 0.064 eV; KATRIN 2025: < 0.45 eV.

**Watch out for.** D3's TT likelihood (83 band powers with independent errors, no foreground or calibration uncertainty, nothing below ℓ = 30) gives errors about 20% too small, and without the low-ℓ data n_s slides and drags N_eff down to 2.6; present it as a demonstration. Planck's plik_lite plus its low-ℓ likelihoods reproduce the published results. DESI DR1's printed Σm_ν < 0.072 eV was corrected to 0.082 after a bug fix in ACT's lensing likelihood. A limit computed with a prior from Σm_ν = 0 cannot be read off the mass plot as a bound on the lightest mass: use a prior starting at each ordering's minimum (normal ordering: lightest < 25 meV, not 6). Cosmological neutrino limits assume a model (ΛCDM) and loosen with evolving dark energy. A "best fit" below zero in some analyses means the data prefer less clustering than ΛCDM allows with any neutrino mass; it does not mean negative mass.

**Builds on.** P2 (the thermal history), D3 (CAMB), D5 (DESI), K1 (BBN and N_eff), C7 (neutrinos).

## P4. Freeze-out: why there is dark matter (if it is a WIMP) ★★☆

**Science question.** If dark matter is a particle that once annihilated in pairs, how many would survive? The answer depends mainly on how strongly they annihilate, and a cross-section typical of the weak force gives about the observed amount: the "WIMP miracle". Can we solve the equation that governs freeze-out, find the cross-section that makes Planck's dark matter density, and see where the searches stand?

**Background.** In equilibrium, the abundance of a heavy particle of mass m falls as exp(−m/T) as the universe cools. When the annihilation rate n⟨σv⟩ drops below the expansion rate H, annihilations stop and the number per comoving volume freezes, at x = m/T ≈ 20–25. The Boltzmann equation `dY/dx = −(s⟨σv⟩ / (H x)) (Y² − Y_eq²)`, with Y = n/s, gives roughly Ω h² ≈ 0.1 × (3 × 10⁻²⁶ cm³/s) / ⟨σv⟩. A precise calculation gives ⟨σv⟩ ≈ 2.2 × 10⁻²⁶ cm³/s for masses above ~10 GeV for WMAP's Ω h² = 0.11 (Steigman, Dasgupta & Beacom 2012); for Planck's Ω_c h² = 0.120 it is 2.1 × 10⁻²⁶ (Saikawa & Shirai 2020). The same physics set the neutron-to-proton ratio for BBN (K1) and the leftover electrons after recombination (O4). Today the same annihilations would make gamma rays in dark-matter-dominated dwarf galaxies. Fermi-LAT sees none and excludes the thermal cross-section for masses below ~100 GeV annihilating into b quarks (Ackermann et al. 2015); with 14 years of data and 42 dwarfs the crossing is ~80 GeV, because of mild excesses (McDaniel et al. 2024). Underground detectors (K6) and the LHC test the same particles other ways. Other candidates are made differently ("freeze-in", or axions).

**Data.** A calculation (SciPy's stiff ODE solvers) with Saikawa & Shirai's g*(T) (`member.ipmu.jp/satoshi.shirai`; P2's own build as a variant), Planck's Ω_c h². Answer keys: Saikawa & Shirai 2020's relic cross-section tables (same page, DM2020/) and Steigman et al.'s ancillary curve (arXiv 1204.3622). Gamma rays: McDaniel et al. 2024's data release on figshare (doi:10.6084/m9.figshare.24058650: per-dwarf, per-energy-bin likelihoods of 14 years of LAT data, their TS profiles and J-factors, 4 MB; the blank-field sets are 4.5 GB each) with PPPC 4 DM ID / CosmiXs photon spectra (gammapy-data mirror on GitHub; the PPPC site refuses scripts); Ackermann et al. 2015's limits read from the vector figures in their arXiv source. Direct detection (J6 and K6 not built): LZ 2024 on HEPData (record ins2841863; needs a non-browser user agent), XENONnT 2025 (GitHub XENONnT/wimp_data_release), C. O'Hare's limit compilation and xenon neutrino floor (GitHub cajohare/NeutrinoFog). For the Higgs-portal relic line, the LHC Higgs WG Yellow Report 4 width table (twiki.cern.ch/twiki/pub/LHCPhysics/HiggsXSBR/Higgs_XSBR_YR4_update.xlsx). For the neutrons, PRIMAT's nuclear rates and the BBN code PRyMordial (GitHub vallima/PRyMordial) as the check; Escudero 2020's background.

**Method + visualizations.**
- Y(x) for several cross-sections: equilibrium, freeze-out and the final plateau (the famous Kolb & Turner plot).
- The ⟨σv⟩ needed for Ω_c h² = 0.120 against mass from 1 GeV to 10 TeV, with the feature where g* changes at the QCD crossover.
- The analytic freeze-out approximation against the numerical solution.
- Neutron freeze-out with the same code: n/p at ~1 s → K1's helium.
- The WIMP status plot: the thermal line against Fermi-LAT's dwarf limits (b b̄ and τ⁺τ⁻) and what mass range remains, with K6's direct-detection plane alongside.
- Extension: stack Fermi-LAT photons on J6's dwarfs (K6's suggested extension). (Done as our own combination of the released per-dwarf, per-bin likelihoods; a photon-level reanalysis needs the Fermitools, which have no Python 3.14 build.)
- Added: the leftover electrons after recombination (Peebles atom vs CAMB), and a complete thermal-relic model on the direct-detection plane (the Higgs-portal scalar singlet).

**Check against.** Saikawa & Shirai 2020: 2.11 × 10⁻²⁶ cm³/s at 100 GeV (Ω h² = 0.120); Steigman et al. 2012: 2.2 × 10⁻²⁶ cm³/s (Ω h² = 0.11); Kolb & Turner (1990)'s analytic freeze-out; McDaniel et al. 2024's TS profiles; Ackermann et al. 2015; PRyMordial (Burns, Tait & Valli 2024) for the neutrons.

**Watch out for.** The equation is stiff at small x; use an implicit solver or start from the equilibrium approximation. The thermal value assumes s-wave annihilation of a self-conjugate particle; Dirac particles, p-wave annihilation, coannihilation and resonances change it. Gamma-ray limits depend on the dwarfs' dark matter content (J-factors), uncertain by factors of a few for the faintest. Dirac conventions differ by 2 between papers (⟨σv⟩ per n_DM², as Saikawa & Shirai, vs the physical χχ̄ cross section). McDaniel et al.'s release notebook evaluates the spectrum at bin centres and holds the likelihood flat beyond the scan; bin-integrated fluxes give 1–2% weaker limits. For the neutron part, the baryon density during e⁺e⁻ annihilation follows the scale factor, not today's η × n_γ (it is 2.7× higher at 1 MeV, 3% when deuterium burns).

**Builds on.** P2 (g*), K6 (direct detection), J6 (dwarfs), K1 (neutron freeze-out), O4 (recombination freeze-out), K5 (the local density).

## P5. Where did the antimatter go? ★☆☆

**Science question.** The laws of physics treat matter and antimatter almost identically, and the Big Bang should have made equal amounts. Yet the universe is made of matter: for every ~270 million quark–antiquark pairs that annihilated, one extra quark survived to make us, and today there is one baryon for every 1.6 billion photons. How do we know there are no antimatter galaxies, how well do two independent clocks of the early universe agree on that small excess, and why does the Standard Model fall short of explaining it?

**Background.** The baryon-to-photon ratio η ≈ 6.1 × 10⁻¹⁰ is measured twice: from the deuterium made in the first minutes (K1) and from the CMB's acoustic peaks at 380,000 years (D3), and the two agree. Antimatter regions are ruled out: where matter and antimatter met, annihilation would make gamma rays, and the diffuse gamma-ray background excludes domains smaller than ~1 Gpc, essentially the whole visible universe (Cohen, De Rújula & Glashow 1998). Colliding galaxy clusters such as the Bullet Cluster (K3) limit antimatter in their gas to a few parts per million (Steigman 2008). Cosmic rays contain antiprotons at about the 10⁻⁴ level expected from collisions (AMS-02), and no confirmed antihelium. Making the excess needs Sakharov's three conditions (1967): baryon-number violation, C and CP violation, and a departure from thermal equilibrium. The Standard Model has all three in principle, but its CP violation (Jarlskog invariant J ≈ 3 × 10⁻⁵) is suppressed by small quark-mass differences to ~10⁻¹⁹–10⁻²³ (depending on the temperature and on which quark masses are used), far below the ~10⁻¹⁰ needed, and its electroweak transition is a smooth crossover (P2) rather than the out-of-equilibrium transition required. Leptogenesis through heavy neutrinos is a leading alternative (P3). Antimatter does fall down: ALPHA-g measured antihydrogen's acceleration at 0.75 ± 0.13 (stat. + syst.) ± 0.16 (simulation) g (Anderson et al. 2023).

**Data.** η from D3 and from deuterium: the PDG's 12 quasar-absorber D/H measurements run through the BBN code PRyMordial (GitHub `vallima/PRyMordial`, not on PyPI), with both its rate compilations; K1 can reuse this. The MeV gamma-ray background from COMPTEL, 0.8–30 MeV (Weidenspointner 1999 and Kappadath 1998 thesis tables), is essential for the domain argument: the glow of annihilations at z > 19 lands at 1–45 MeV, below Fermi-LAT's band. Fermi-LAT's isotropic gamma-ray background (Ackermann et al. 2015, VizieR J/ApJ/799/86, foreground models A–C) and EGRET (Strong et al. 2004) cover higher energies. From the Cosmic-Ray DataBase (CRDB, `lpsc.in2p3.fr/crdb`, which needs a non-browser User-Agent): AMS-02's antiproton-to-proton ratio, plus B/C and the p, He, C, N, O, Ne, Mg and Si fluxes for a leaky-box grammage. p̄p annihilation at rest uses Crystal Barrel channel fractions (Klempt, Batty & Richard 2005). CKM parameters come from the PDG. Published bounds: Cohen et al. 1998, whose figures are vector PostScript in their arXiv source and readable exactly, and Steigman 2008 with the Fermi-LAT cluster limits of Ackermann et al. 2010.

**Method + visualizations.**
- η two ways: BBN deuterium (K1) and the CMB (D3) on one axis, converted between η and Ω_b h².
- The asymmetry in numbers: today's 411 photons per cm³ against ~0.25 baryons per m³.
- The domain argument: the gamma rays expected from annihilation at domain boundaries (the pion-decay bump at ~70 MeV, redshifted) against COMPTEL's MeV background (Fermi-LAT and EGRET above 30 MeV) → a lower limit on domain size.
- Cosmic-ray antiprotons: AMS-02's p̄/p against energy, beside a simple estimate of the antiprotons made by cosmic-ray collisions.
- The Sakharov conditions and why the Standard Model fails: the CP-violation measure against what is needed, and P2's crossover.

**Check against.** Planck 2018 and Cooke et al. 2018 (K1): η ≈ 6.1 × 10⁻¹⁰; Pitrou et al. 2021 (PRIMAT rates): Ω_b h² = 0.02195 ± 0.00022 from Cooke's D/H; Cohen, De Rújula & Glashow 1998: domains ≳ 1 Gpc. That bound is conservative: their own curve against COMPTEL's final data needs ~10 Gpc, and today's cosmology raises it further. AMS-02 2016: p̄/p ≈ 2 × 10⁻⁴ at high energy.

**Watch out for.** The domain bound assumes the domains were in contact after recombination; be explicit about the assumption. It also depends strongly on when annihilation stops (structure formation) and on how much annihilation electrons heat the gas, so show both. Antiproton predictions depend on cosmic-ray propagation models, and claims of an excess around 10–20 GeV are contested. Korsmeier et al. 2018's eq. 14 prints (10 − √s)² in the low-energy factor of their Param. II; only the first power (Winkler 2017) reproduces their own Fig. 4, and the square fakes a 2.5× p̄/p bump at ~9 GV. The two BBN rate compilations differ by ~2% in D/H, which shifts η by ~1.3σ. Baryogenesis models are mostly unconstrained by data: this project shows what the data require, not which model is right.

**Builds on.** K1 (BBN), D3 (CMB), K3 (the Bullet Cluster), P2 (the electroweak crossover), P3 (neutrinos).

## P6. Recreating the first microseconds: the quark–gluon plasma ★★★

**Science question.** About 10 microseconds after the Big Bang, the universe was a plasma of free quarks and gluons at nearly 2 trillion kelvin. The LHC and RHIC recreate it by smashing lead or gold nuclei together. Can we take the particle yields measured by ALICE, fit them with a thermal model, and find the temperature at which quarks bound into hadrons, the same transition the whole universe went through?

**Background.** A central lead–lead collision at the LHC sends ~1,600 charged particles into each unit of rapidity. Their relative abundances (pions, kaons, protons, Λ, Ξ, Ω, φ, deuterons, helium-3 and their antiparticles) match a hadron gas in chemical equilibrium at a single temperature, once decays of heavier resonances are added. Each primary yield scales roughly as (2J + 1) m² T K₂(m/T) × e^(μ/T) (Boltzmann limit). The fit to ALICE data gives T = 156.5 ± 1.5 MeV, with a baryon chemical potential near zero, as in the early universe (Andronic et al. 2018). That matches lattice QCD's crossover temperature, 156.5 ± 1.5 MeV (HotQCD 2019), which the universe passed at ~10–20 µs (P2). The plasma starts hotter: direct photons measured by ALICE point to ~300 MeV (T_eff = 297 ± 12 (stat) ± 41 (syst) MeV after subtracting prompt photons, 304 ± 11 ± 40 without; ALICE 2016).

**Data.** ALICE Pb–Pb yields at √s_NN = 2.76 TeV on HEPData: π, K and p (INSPIRE 1222333) plus the published K⁰_S, Λ, Ξ, Ω, φ, deuteron and ³He measurements (HEPData ins1243863, ins1243865, ins1288320, ins1380491, ins1380234; the deuteron, ³He and ⁴He yields themselves are only in the papers' tables, in the arXiv TeX sources; the hypertriton is published as yield × branching ratio, which ALICE took as 25%; Λ̄ was never published separately). Thermal-FIST's file `input/data/ALICE-PbPb2.76TeV-0-10-all.dat` holds exactly the values Andronic et al. fitted. STAR's beam-energy-scan yields (HEPData ins1510593, ins1738953) allow the same fit at 7.7–39 GeV. HEPData's Cloudflare front refuses browser and Python-default user agents; a plain `python-requests` one works. Hadron masses, spins and decay branching ratios from the PDG via the `particle` package, or the particle and decay lists shipped with the open-source Thermal-FIST code (GitHub `vlvovch/Thermal-FIST`). Lattice QCD results (HotQCD 2019; Borsányi et al. 2020).

**Method + visualizations.**
- A hadron resonance gas: every hadron up to ~2 GeV as an ideal quantum gas → primary densities at temperature T and chemical potentials.
- Feed-down: decay each resonance into stable daughters using the PDG branching ratios → final yields.
- Fit T, μ_B and volume to central-collision yields → the famous plot of measured against predicted yields across nine orders of magnitude, from pions to helium-3 and its antinucleus.
- The QCD phase diagram: T against μ_B with ALICE's point, lower-energy experiments' freeze-out points (published fits), the lattice crossover band, and the early universe's path along μ_B ≈ 0.
- Lattice QCD's equation of state (pressure and energy density over T⁴) against the hadron gas below T_pc and a free quark–gluon gas above it: the transition seen in P2's g*.

**Check against.** Andronic et al. 2018 (*Nature*): T = 156.5 ± 1.5 MeV, μ_B = 0.7 ± 3.8 MeV; HotQCD 2019: T_pc = 156.5 ± 1.5 MeV.

**Watch out for.** Feed-down dominates the pion and proton yields, and experiments differ in how they count weak-decay products (Λ → pπ); follow ALICE's corrections. The proton yield sits below the fit (the "proton anomaly"), and short-lived resonances such as K* are depleted by rescattering after freeze-out; exclude or discuss them. Loosely bound nuclei matching a 156 MeV source is a puzzle in itself. Andronic et al.'s volume (5,280 fm³) assumes an excluded volume of r₀ = 0.3 fm for every hadron; with point-like hadrons it is 24% smaller, and T and μ_B are unchanged. They also add the pions and kaons from charm-hadron decays, which ALICE counts as primary; this raises the fitted T by ~0.9 MeV. Model implementations agree to < 1% for species without feed-down but differ by 3–15% where many resonances decay into a species (Λ, φ), worth ~0.6 MeV in T.

**Builds on.** P2 (the QCD epoch), C8 (dense matter in neutron stars, the other corner of the phase diagram), L4 (nuclei).

## P7. The Higgs field switches on: Higgs bosons in CMS open data ★★☆

**Science question.** About 10⁻¹¹ seconds after the Big Bang the Higgs field settled into its present value, giving W and Z bosons, quarks and electrons their masses. Its particle, the Higgs boson, was discovered at CERN in 2012. Can we find it ourselves in CMS's public 2012 data, see the particle zoo as bumps in a mass spectrum, and ask what the measured Higgs and top-quark masses say about the stability of the vacuum we live in?

**Background.** A Higgs boson occasionally decays into two Z bosons, each of which decays into an electron or muon pair. The four leptons' combined invariant mass, `m = √((ΣE)² − |Σp|²c²) / c²`, peaks at the Higgs mass, 125.13 ± 0.11 GeV (PDG 2026; 125.20 in PDG 2024), above a smooth background of directly produced Z pairs. This channel and the two-photon one carried the 2012 discovery (ATLAS and CMS). The same method applied to muon pairs reveals the particle zoo, from the η and ρ/ω mesons through the J/ψ and ϒ to the Z at 91 GeV; each was a species the early universe made and lost as it cooled (P2). The Higgs self-coupling λ changes with energy (renormalization-group running), dragged down by the heavy top quark. With the measured masses (m_t = 172.60 ± 0.27 GeV, PDG 2026; 172.57 ± 0.29 in PDG 2024) λ turns slightly negative around 10¹⁰–10¹¹ GeV: our vacuum is probably metastable, with a lifetime vastly longer than the age of the universe (Degrassi et al. 2012; Buttazzo et al. 2013).

**Data.** The CERN Open Data portal's CMS 2012 data in reduced NanoAOD format for education: the Higgs-to-four-leptons analysis (record 12360) and its datasets: simulated Higgs (12361, 42 MB), simulated ZZ background (12362–12364, ~0.5 GB), and the DoubleMuParked and DoubleElectron data (12365–12368, ~12 GB in all). The dimuon-spectrum analysis (12342) and its dataset (12341, 2.2 GB). Read with `uproot` and `awkward` (wheels exist for Python 3.14); the files also serve byte ranges over HTTPS (`opendata.cern.ch/eos/opendata/...`). Streaming only the needed branches saves only ~20%, so it is simpler to fetch whole files (15.0 GB in all, ~1.2 MB/s per connection, so use parallel byte ranges) and check them against the adler32 checksums the portal lists. The outreach code's published plots are vector PDFs (in the record's zip) and can be read back bin by bin. PDG masses; two- and three-loop Standard Model β functions and NNLO matching at the top mass from Buttazzo et al. 2013 (appendices A and B, in the arXiv source); SMDR (Martin & Robertson 2019, C, builds with `-std=gnu17`) as an independent check.

**Method + visualizations.**
- The dimuon spectrum on a log–log scale with every resonance labelled: the particle zoo in one plot.
- The four-lepton analysis: lepton selection, Z-pair reconstruction, and the m₄ℓ histogram with the Z → 4ℓ peak at 91 GeV and the Higgs at 125 GeV over the simulated ZZ background → our own Higgs mass and an indicative significance.
- The Higgs potential: the "Mexican hat", the 246 GeV vacuum value, and the W, Z and fermion masses it gives.
- Running couplings: integrate the Standard Model's renormalization-group equations for the gauge couplings, the top Yukawa coupling and λ from the measured masses up to the Planck scale → where λ crosses zero; the stability diagram (m_t against m_H) with the measured point in the metastable band.
- On P2's calendar: when the electroweak crossover happened and what it changed.

**Check against.** The outreach analysis's published plots (a visible 125 GeV excess in this dataset), bin by bin; ROOT's df102 tutorial cut flow for the dimuon file (61,540,413 → 31,104,343 → 24,067,843); PDG 2026 masses; Buttazzo et al. 2013 (their table of couplings at M_t, the couplings at the Planck mass, instability scale ~10¹⁰–10¹¹ GeV, metastable vacuum); Andreassen, Frost & Schwartz 2018 for the decay rate.

**Watch out for.** This open-data analysis has only a handful of Higgs candidates and ignores the minor backgrounds, so any significance is indicative; the official CMS analysis used more data and event categories. The example code applies no lepton efficiencies, only a global normalization (luminosity × σ/N, ×1.386 for ZZ): check it in the on-shell ZZ region above 180 GeV. Its 2e2μ data add the DoubleMuParked and DoubleElectron streams, which overlap (54 of 94 collisions, three of them in the Higgs peak): remove duplicates by run, luminosity block and event. The files' lepton momenta are uncalibrated (muons ~0.15% low, electrons 0.4–1% high relative to the simulation): calibrate them on the Z peaks before measuring m_H. One-loop running puts the instability at the wrong scale; use two loops, and show how the answer moves with m_t, which dominates. "Metastable" is a statement about the Standard Model with no new physics up to enormous energies.

**Builds on.** P2 (the electroweak epoch), P5 (the electroweak crossover and baryogenesis), C7 (particle physics in astronomy).

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
| A4 Earth as an exoplanet | **done** — [44-earth-as-exoplanet/](44-earth-as-exoplanet/). All 55,179 public EPIC natural-colour images (June 2015 – October 2026) turned into one point of light each; calibrated level-1B data need a NASA Earthdata Login and were not used (spec deviation, documented). Day length: 23.993 ± 0.004 h (2015–16) and 24.008 ± 0.003 h (2024) from one year; **24 h − 0.15 ± 0.28 s from eleven years**, the solar day (sidereal excluded: power 0.006 vs 0.29), since DSCOVR orbits the Sun with Earth. Jiang et al. 2018's sampling test redone: ≥90% detection to ~5 h between looks at random times, aliasing below ~4 looks a day on a fixed schedule. The dot's colour tracks the land share of its light (r = −0.95); the Pacific side is ~17% bluer and ~12% fainter in red than the Africa side; colour variations are 96% geography, brightness ~40–60% weather; PC2 = surface (r = 0.96 with land), as Fan et al. 2019. Seasons as Yang et al. 2018 (brightest December, second peak May–June). Phase angle 1.8–12.4° (not 4–15°) on an 87.5-day cycle; brightening ~0.9%/° toward small angles, ~4× a Lambert sphere. Moon: 16 transits and 22 occultations predicted from JPL Horizons (DSCOVR = body −78); EPIC imaged 31 of 35 into a login-gated collection, three transits are public. The display images' unpublished tone curve (six tests of gamma, frame-to-frame jumps of 1.5% rms) limits every amplitude to a factor ~1.5 and the Moon/Earth ratio to ~2. Plus the keeper: a year of 2024 as a movie with its dot and light curve (30 MB), and transit GIFs. Spec corrections above: day count, Earthdata Login, phase-angle range, 2025 gap. |
| A5 ALMA protoplanetary discs | **done** — [45-protoplanetary-discs/](45-protoplanetary-discs/). The 20 DSHARP discs and HL Tau measured from the released images by our own pipeline (2.5 GB). Orientations from MCMC ellipse fits to traced ring crests, or beam-deconvolved Gaussians: all 18 measurable inclinations within 2σ of DSHARP II (AS 209 35.05° vs 34.97°; HL Tau 46.85° vs ALMA 2015's 46.72°; nearly face-on Sz 114 unconstrained). Profiles match DSHARP's released ones to < 1% for 13 of 18 discs; 29 of 34 ellipse-traced rings and gaps recovered at a median 1.6 mas (0.04 beam), widths within 0.13 au, depths within 1%; DSHARP's "plateau" features (16) are not extrema and are not recovered. DSHARP VII's width-to-mass relation on our widths reproduces its planets (median ratio 0.99, 11 of 12 gaps): 0.03–2.2 M_J at 9–99 au, plus 19 M_J in HD 143006's inner gap. Parsing the papers' tables showed DSHARP VII's Elias 24 masses imply Δ = 0.334, not the printed 0.32, and its AS 209 "99 au" gap is D90 + D105 merged; MY Lup's released profile uses 2-au bins. Gap temperatures show no snow-line clustering. HD 163296 weighed with its CO velocity field (flared Keplerian surface; all four sky orientations tested, near side NE): 2.058 ± 0.034 M☉ at 101 pc (Teague et al. 2019: 2.022), emission surface z₀ = 284 vs 265 mas; the rotation curve's pressure minima at 81 and 139 au match Teague et al. 2018's two planets (83, 137 au); Pinte et al.'s kinks are not seen in line-centre residuals (0.8–0.9σ). Spiral pitch angles by arm tracing and a log-spiral Fourier transform (WaOph 6 13.9°/14.2° vs 13.4°/13.7°). Spec corrections above: download 2.5 GB, not 1.3 GB; HL Tau and Elias 27 need cloud distances. Plus an offline unroll-a-disc explorer. |
| A6 Directly imaged orbits | **done** — [46-hr8799-orbits/](46-hr8799-orbits/). 167 positions of HR 8799 b–e from 1998–2025 (HST, Keck, GPI, SPHERE, LBT, GRAVITY, JWST), parsed from the papers' LaTeX and fitted jointly by our own nested sampling (sky-state parametrisation, Jacobi geometry, per-instrument north/scale offsets and jitter). Wang et al. 2018 reproduced: every element of their four models within 1.43σ; the Bayesian evidences they could not compute favour near-1:2:4:8 coplanar orbits by ln B = 18.6 (25.9 with all data). The planets re-derive GPI's north-angle error (+0.34 ± 0.09° vs De Rosa et al.'s +0.36°). Orbits fitted to 2004–2016 predict 77 later positions to 9 mas (GPI recalibrated); the final orbits predict 17 unpublished GRAVITY 2019–2022 positions found in whereistheplanet's data file to 2.5 mas. N-body (REBOUND, Wang's set-up, 470,000 draws): 191 of 200,000 stable for 40 Myr (Wang: 441 per million), matching their Table 3 within 0.5σ, M_cde < 7.3 MJup; hot-start masses 1 of 150,000; the all-data orbits are stable 16× more often (1.6%; 1% even without the period-ratio constraint), and their survivors land on whereistheplanet's GRAVITY-informed orbits. Resonances as Wang found (d:e locked, Laplace rare, b free); chaos: a 15-cm nudge of b, or a last-bit change of a mass (a pandas CSV round-trip bug, found and fixed), breaks up half the survivors; 34% last another 100 Myr. HGCA: m_e = 9.9 ± 2.0 MJup (Brandt et al.: 9.6). β Pic b: closest approach 2017-09-12 ± 3 d at 27 R★, Hill sphere in front Apr 2017–Feb 2018 (Nielsen et al. 2020). orbitize!: OFTI fails for d and e, its MCMC got trapped for β Pic. Keeper animation 1998–2050 and an offline orbit explorer. Spec corrections above: download ~115 MB (1 MB kept), periods ~500/55 yr, orbitize! as the cross-check. |
| B1 HR diagram | **done** — [02-hr-diagram/](02-hr-diagram/). 215,844 Gaia DR3 stars within 100 pc. Twin-binary sequence 0.75 ± 0.02 mag above the main sequence (0.753 predicted); red clump M_G 0.51 ± 0.02; the Jao gap confirmed at 5.6σ; the white dwarf A/B split; Parenago's discontinuity at BP−RP 0.84. Plus an all-sky star-count map of 1.8 billion sources. |
| B2 3D solar neighborhood | **done** — [06-solar-neighborhood/](06-solar-neighborhood/). 300,338 stars within 100 pc: 299,525 from Gaia DR3 (72,701 spurious sources removed with the Rybizki et al. fidelity classifier, 71% of them toward the crowded inner Galaxy) plus 813 from Hipparcos that are too bright or too binary for Gaia. Interactive 3D map in true star colours, with SIMBAD names for every star within 20 pc. Red dwarfs are 74% of the stars within 20 pc and none of the 205 visible to the naked eye; only 46% of the naked-eye sky lies inside 100 pc. The Hyades found blind as the densest clump, at 46.8 pc (literature 47.5), and Coma Berenices at 85.5 pc (85.9); both beat an isotropic null (142 and 48 neighbours vs 25). Photometric classes match SIMBAD spectral letters for 96% of dwarfs; 99.3% of the Gaia Catalogue of Nearby Stars recovered. Of 477 objects known within 10 pc, 330 found: every missed star is in a binary, and 70 of the 79 missed brown dwarfs are too faint. 2.0 ± 0.5% more stars below the Sun than above (4.4σ), the sign expected for a Sun ~20 pc above the disc. |
| B3 Cluster hunting | **done** — [02-cluster-hunting/](02-cluster-hunting/). The Pleiades, Praesepe and M67 found blind as knots in Gaia DR3 proper motions (1,358 / 953 / 1,208 members); distances within 0.2% of Hunt & Reffert (2023), Pleiades at 135.0 ± 0.2 pc; members' main sequences 2.3–6.5× thinner than the field's; perspective and depth effects of a shared 3D velocity confirmed. Isochrone age fitting is the open extension. |
| B4 Cepheids / Leavitt law | **done** — [47-cepheids/](47-cepheids/). All 9,667 OGLE-IV Magellanic Cepheids, periods and intensity means measured blind from the raw light curves (OGLE's period recovered for 99.9%, means to 0.0004 mag). LMC Wesenheit law −3.316 log P + 15.890, scatter 0.078 mag (0.207 in V); all twelve OGLE fits (Soszyński et al. 2015) reproduced. Leavitt's 1912 table transcribed from the scan: her slope −2.03 ± 0.10 (Breuval et al. 2025: same); all 25 stars found in OGLE, three of her periods aliases of one cycle per year or per lunar month. LMC tilted disk i = 25.8° ± 1.0°, line of nodes 148° ± 2°, north-east nearer; SMC 3.2 kpc deep, SMC − LMC 0.537 ± 0.006 mag. Zero point from 630 Milky Way Cepheids with Gaia DR3 parallaxes and spectroscopic [Fe/H] in a parallax-space fit with the LMC distance as a parameter: **μ_LMC = 18.495 ± 0.041 ± 0.053 (50.0 ± 1.5 kpc)** vs eclipsing binaries 18.477 ± 0.026; metallicity term −0.32 ± 0.07 (SMC − LMC test in Gaia bands: −0.26 ± 0.05), residual parallax offset 22 ± 3 µas, dust coefficient 1.76 (not 1.90), Gaia parallax errors ×1.8. Found and corrected a 0.036 mag regression-dilution bias in the standard method (100 injection simulations). SH0ES with our LMC alone: H0 = 73.0 ± 2.2; Planck would need the LMC at 54 kpc (2.6σ). Plus a Petersen-diagram metallicity sequence and an offline 3D map of 9,365 Cepheids. |
| B5 Asteroseismology | **done** — [48-asteroseismology/](48-asteroseismology/). Our own pipeline (Parseval-calibrated spectra, blind maximum-likelihood background + envelope fit, autocorrelation Δν refined with the universal pattern) on 4 years of Kepler photometry of 104 red giants (νmax 6–211 µHz): νmax +0.23% and Δν −0.12% from APOKASC-3 with 1.0% and 0.6% scatter, mid-pack among their ten pipelines; masses on APOKASC-3's scale +1.5% (scatter 3.3%), radii +0.5% (1.4%). The Sun through the same code (SOHO/VIRGO, 26 years): νmax☉ 3142.6 ± 3.1 µHz (SYD uses 3090), rising 27 ± 7 µHz with activity (Howe et al. 2020: ~25). Uncorrected radii vs Gaia: +1.9% (RGB) / −1.8% (clump) with our Sun, +4.8% on the RGB with 3090. Own stretched-period ΔΠ₁ for 41 giants: shell- vs core-burning classification 41/41 with APOKASC-3, 28 of 29 values within 2% of Vrard et al. 2016. 16 Cyg A/B radii within 1–2% of CHARA interferometry. Sound files of the real light curves and an interactive échelle page. Spec's download estimate corrected (~1.1 GB incl. 0.23 GB of VIRGO). |
| B6 Oldest stars | **done** — [49-oldest-stars/](49-oldest-stars/). HUGS photometry of M92, NGC 6397, M15 and 47 Tuc fitted with three α-enhanced isochrone sets (MIST v2.5, BaSTI-IAC, DSED). Distances from Baumgardt & Vasiliev (2021) with a 0.025 mag systematic floor; 1,000-draw Monte Carlo over distance, reddening, [Fe/H] and [α/Fe] per set. Ages: M92 13.42 ± 0.76 Gyr, NGC 6397 13.98 ± 1.01, M15 12.84 ± 0.72, 47 Tuc 11.36 ± 0.90, which formed 2.0 Gyr after the metal-poor three. M92 matches Ying et al. 2023 (13.80 ± 0.75) once their 0.07 mag larger distance modulus is used (Dartmouth models: 13.4–13.8). No single error term dominates: distance, reddening, composition and model set are each 0.2–0.6 Gyr. DSED runs 0.6–1 Gyr older than MIST/BaSTI for the same diagram. The oldest clusters plus formation time need a universe ≥ 13.6 ± 0.8 Gyr old: consistent with D3 (13.83) and D8 (13.45), 1.3σ above D1's 12.46 (H0 = 73.5); H0 < 75.5 km/s/Mpc (95%) for D3's Ωm, < 48 without dark energy (the 1990s age crisis). Chromosome maps split every cluster into two generations; first-generation fractions agree with Milone et al. 2017 to ~0.1. Extension, all 56 HUGS clusters (29 clean): ages 0.6–0.7 Gyr older than VandenBerg 2013 and Valcin 2025; in-situ clusters 1.1 Gyr older than accreted ones at [Fe/H] −1.7 to −1.0; fitted colour shifts grow with reddening (4σ), so either MIST's E(F606W−F814W)/E(B−V) = 1.17 should be ~1.0–1.1 or the Harris reddenings are ~10% high. Traps found: the recommended HUGS quality cut deletes every single-exposure (bright) star, whose RMS is 0; MIST v2.5's ACS link is broken. Offline fit-it-yourself page. Download ~2.1 GB, +2.6 GB for the survey (spec said ~100–300 MB). |
| B7 Spectroscopy | not started |
| B8 Cluster time-lapse | not started |
| B9 Eclipsing binaries | not started |
| B10 White dwarfs | not started |
| B11 MESA to core collapse | not started |
| B12 Constellations in motion | not started |
| B13 Stellar atmospheres | not started |
| B14 Starspots and activity cycles | not started |
| B15 Stellar coronae (X-rays) | not started |
| B16 Stellar winds | not started |
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
| D5 Evolving dark energy | **done** — [42-evolving-dark-energy/](42-evolving-dark-energy/). Our own likelihood (DESI DR2 BAO, four supernova samples, compressed Planck PR4 prior; expansion history checked against CAMB to 10⁻⁶) reproduces DESI's released best-fit χ² to 10⁻⁴, its Δχ² to 0.02 and its chains to 0.08σ. ΛCDM: Ωm = 0.2978 ± 0.0086 (DESI 0.2975). w0wa over ΛCDM: DESI 1.7σ, DESI + CMB 2.4σ; with supernovae 2.3σ (Pantheon+), 3.3σ (Union3), 3.8σ (DESY5) and 2.8σ (DES-Dovekie, new). The evidence tracks how far each sample's ΛCDM Ωm sits from DESI's (1.7σ–2.9σ). Found: DESI's compressed CMB prior compares ω_b + ω_c with a mean that includes the neutrino density; a consistent mean lowers those results by 0.1–0.3σ (the full-CMB headline numbers are unaffected). Our own BAO in 1.28 M DESI DR1 LRGs (own C pair counter and FFT reconstruction): peak at 6.6σ / 4.5σ after reconstruction (4.2σ / 3.1σ before), D_V/r_d = 12.44 ± 0.18 (z = 0.51) and 15.86 ± 0.21 (z = 0.70), plus D_M/r_d and D_H/r_d separately, all within 1σ of DESI DR1. Plus an interactive w0–wa profile-likelihood explorer. |
| D6 Fresh supernovae | **done** — [39-fresh-supernovae/](39-fresh-supernovae/). All 3,628 ZTF SN Ia DR2 light curves (only 4 in Pantheon+) fitted twice with sncosmo: SALT3 (matches the ZTF team's own SALT3 fits to 0.034 in x1, 0.005 in c, 0.003 mag in mB) and D1's own SALT2 model from the Pantheon+ release. D1's α = 0.148, β = 3.09 are recovered exactly from D1's table (with an emulated Pantheon+ bias correction). The trap: SNANA adds the model's `MAG_OFFSET` of 0.27 mag and sncosmo does not (measured 0.262 ± 0.002 on 40 of Pantheon+'s own Foundation light curves); ignoring it puts every supernova 0.19 mag too faint. Done D1's way, 2,389 fresh supernovae land −0.027 ± 0.006 mag from D1's curve (0.023 < z < 0.06), well inside ZTF's blinded, ~5%-calibrated flux scale, so this is not an H₀ measurement. Scatter nMAD 0.151 mag (ZTF team 0.150; 0.318 before standardisation), with no leftover stretch trend or host-mass step (bare Tripp: 0.087 mag step). ZTF-fitted β = 3.07 ± 0.04 (Ginolin+25a 3.05). Selection bias beyond z = 0.06 reproduced by a forward model from the published selection function; spectral-feature redshifts need σz = 0.0037; random galaxy velocities not detected (< 230 km/s). Open: the nearest ZTF supernovae (z < 0.03) are ~0.08 mag faint; a dipole fits 444 ± 93 km/s towards (304°, 13°), Carrick+15's bulk-flow direction, but they are faint on both sides of the sky, so a flux non-linearity is not excluded. Live: 25 SNe Ia classified on TNS 20 Aug–4 Oct 2026 (ZTF alert forced photometry via ALeRCE) land at −0.06 ± 0.04 mag. Offline explorer to standardise them yourself. |
| D7 Is the universe expanding? | **done** — [37-expansion-tests/](37-expansion-tests/). Three tests; tired light fails each. **Time dilation** (DES-SN5YR): White et al.'s (2024) reference-curve method replicated on the same 1,504 supernovae, b = 1.002 ± 0.006 (published 1.003 ± 0.005), per-supernova widths matching theirs to ~1%; their Nelder–Mead start at the expected answer makes no net difference. Our own global template fit on the raw photometry (SALT3's fitted points sit in a rest-frame window that already assumes (1 + z), so ours uses a fixed observer-frame window) gives b = 1.012 ± 0.025 (bootstrap; 3.3× the χ² curvature error, from per-supernova stretch), Δχ² = 14,099 against b = 0, per-redshift-bin factors that follow (1 + z), and recovers injected b = 0, 0.5 and 1. **CMB temperature**: 46 published measurements to z = 6.34 give β = 0.004 ± 0.011 (tired light χ² = 5,346 against 23); T₀ from the distant data alone 2.7232 ± 0.0091 K; Noterdaeme+11 (−0.005 vs −0.007) and Riechers+22 (0.0032 vs 0.0034) reproduced; all 107 transcribed values checked against the papers, finding a typo in Klimenko+20's summary table. **Tolman test**: Lubin & Sandage (2001) reproduced to 0.04 mag per cluster (n = 2.59 R, 3.37 I); flat ΛCDM gives 2.29 ± 0.29 and 3.05 ± 0.21, needing the expected brightening of younger stars; with tired-light radii 1.44 and 2.19, which would need galaxies fainter in the past (1.6σ and 5.7σ before any evolution; Paper IV's tired-light case kept expansion's local (1 + z)⁴ correction). |
| D8 Cosmic chronometers | **done** — [38-cosmic-chronometers/](38-cosmic-chronometers/). The 32 cosmic-chronometer H(z) measurements (Moresco et al. 2022, 0.07 < z < 1.97) with their full stellar-population covariance. Published fits reproduced (primer notebook H0 65.99 vs 66.0; review 66.7 ± 5.4 vs 66.5 ± 5.4, Ωm 0.334 vs 0.34), which requires applying the covariance to all 32 points. Found that the published recipe scales the shared systematic by each measured H value, which biases fits low (Peelle's puzzle): −5.7 km/s/Mpc on average in 300 mock data sets, and adding the stellar-library term moves the real H0 from 66.7 to 62.2. Scaling by a smooth reconstruction removes it. Chronometers alone, flat ΛCDM: H0 = 69.0 ± 5.4, Ωm = 0.319 +0.072 −0.058; the expansion turned from slowing to speeding up at z_t = 0.62 ± 0.16 (D1: 0.59 ± 0.04), q0 = −0.52, age 13.4 +0.7 −0.6 Gyr. Gaussian-process reconstruction: dip at z = 0.59 ± 0.17, accelerating today in 92% of curves. Method checks on exact ΛCDM curves: the GP recovers z_t but shrinks H0 by ~3 and q0 toward zero, and Moresco et al. 2016's two-line method returns z_t = 0.34 for a true 0.61, which explains their low 0.4 ± 0.1. DESI DR2 D_H/r_d with D3's sound horizon shows the dip directly (χ² = 4.0/5 against the chronometers). Moresco 2012's SDSS points re-read from their D4000 medians (74.7, 75.9 vs 74.9, 75.0), and simple cuts on J3's SDSS galaxies give 18–154 km/s/Mpc: selection is most of the method. Nine 2022–2026 measurements all within ±0.9σ. Futures to +100 Gyr: ΛCDM grows 290×; DESI's w0–wa extrapolated stops accelerating in 9.6 Gyr; w = −1.24 (D1's 95% limit) rips in 48.7 Gyr, w = −1.5 in 23.1 Gyr, with Caldwell et al.'s timeline of destruction recomputed. |
| D9 Lensed quasar time delays | **done** — [43-lensed-quasar-delays/](43-lensed-quasar-delays/). 13 delays in six lensed quasars from the public COSMOGRAIL light curves with our own estimator (one intrinsic spline shifted per image + a microlensing spline per image, delays by Bayesian evidence because plain χ² rewards pushing an image into seasonal gaps; knot step chosen by BIC; errors from 3,600 simulated light curves plus the spread over six settings). All 13 agree with the published delays within 1σ: RX J1131 A→D 90.9 ± 2.6 d (Millon+2020: 92.5), HE 0435 A→D 14.1 ± 1.3 (13.8), WFI 2033 B→C 57.5 ± 2.4 (59.4), SDSS J1206 113.8 ± 3.8 (111.8), DES J0408 A→B 114.3 ± 4.7 (112.1). PyCS3 run on the same curves agrees to 1–3 d; I2's ICCF works for short delays but is 10 d short for RX J1131's D. H0LiCOW reproduced from its released posteriors: 73.29 (−1.84/+1.68) (Wong+2020: 73.3 −1.8/+1.7), every lens within 0.1. With our delays swapped into the published lens models: H0 = 73.6 (−1.9/+1.7); PG 1115 alone moves (81 → 96 ± 13, one-season B delay 2 d shorter). Mass-sheet degeneracy shown with our own elliptical power-law fits to HE 0435's HST positions (slopes 1.7–2.3 all fit, H0 = 47–91) and TDCOSMO 2025's chains (H0 = λ_int × 73.9, corr 0.82); our delays change TDCOSMO 2025's 71.6 (+3.9/−3.3) by a factor 0.9999. Interactive delay finder page. Download ~190 MB (spec said <50 MB): the H0LiCOW distance chains. |
| D10 ISW effect | **done** — [40-isw-effect/](40-isw-effect/). Planck SMICA (D3's maps) against three superstructure catalogues, BOSS DR12 and NVSS, every analysis choice frozen with a hash before the real sky was stacked; errors from 1,000 simulated skies. NVSS × CMB: A = 1.01 ± 0.40 (ΛCDM = 1; Planck 2015: 0.95 ± 0.36); all cross-spectra 0.94 ± 0.39 (2.4σ from no dark energy). Granett et al.'s 2008 stack replicates (9.0 ± 2.4 µK; theirs 9.6 ± 2.2) but drops to 2.65σ once the search over radius and N is priced on noise skies, and is not ISW: our own ΛCDM prediction from a Wiener-filtered 3D map of 1.24 million BOSS galaxies is +0.06 µK along those sightlines (the supervoids sit in slightly overdense surroundings), and a matched filter gives A = 0.0 ± 1.4. Hotchkiss et al.'s DR7 stacks reproduced within 0.4σ; 1,228 BOSS voids A = −0.07 ± 0.67; BOSS ISW template −0.03 ± 0.44. Cold Spot −4.62σ in a 5° wavelet (p = 0.006, a posteriori; 16% with a Gaussian filter); the Eridanus supervoid's ISW (−20 µK) explains ~15% of its depth and 3% of its wavelet signal. Traps found: BOSS's sparse early LOWZ chunks, one-pixel cracks when rasterising mangle footprints at Nside 1024, the 1/k² potential of shot noise, and rotated-map errors that halve near the Galactic pole. Plus an offline look-elsewhere playground. |
| D11 Cluster counts | **done** — [41-cluster-counts/](41-cluster-counts/). eRASS1's 5,259-cluster cosmology sample (Bulbul et al. 2024) fitted with our own forward model. The selection function and footprint are not public: rebuilt from Ghirardini et al.'s Fig. 4 curves (extracted from the vector PDF), Clerc et al.'s 2000 s curves, the exposure-map HiPS and Bulbul et al.'s Fig. 2 (13,114 vs 13,116 deg²). At the published best fit the model predicts 4,978 clusters; 4,971 observed (purity-weighted). With the published count rate–mass relation free within its stated widths: S8 = 0.885 (−0.033/+0.037), Ωm = 0.282, σ8 = 0.914 (published 0.86 ± 0.01, 0.29, 0.88); relation fixed: S8 = 0.866 ± 0.003. The clusters side with the CMB (Planck 0.832; D3 grown to today with CAMB 0.845 ± 0.029), not DES Y6 (0.789). Read with a 24% lower (hydrostatic-like) mass scale, the same counts give S8 = 0.78. Without dark energy σ8 would be 1.03, with 3.5× more halos above 10¹⁵ M☉. Per-cluster Planck/eRASS1 mass ratios range 0.35–1.26 by estimator (factor-2.7 scatter in the relation). |
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
| F5 Sunspot anatomy | not started |
| F6 Magnetic cycle and polar reversals | not started |
| F7 Flares and superflares | not started |
| F8 Solar wind | not started |
| F9 CMEs and space weather | not started |
| F10 Million-degree corona | not started |
| F11 Helioseismology | not started |
| F12 Solar irradiance | not started |
| G1 Galaxy composites | **done** — [13-galaxy-composites/](13-galaxy-composites/). Twelve galaxies (M51, M101, M74, M63, M81, M64, M104, NGC 5866, NGC 4631, M82, the Mice, M87) built from calibrated SDSS g, r, i with one physical recipe: Lupton et al.'s arcsinh stretch fixed in surface-brightness units, and white = the Sun's measured colour (the survey viewer's default whitens a colour no star has). Checked first: the survey's own JPEGs reproduced (99.995% of 7.2 M values within ±1); 278 stars within 0.004 mag of Gaia's synthetic photometry; whole galaxies within 0.006 mag of Brown et al. 2014 (colours 0.003) once foreground stars are filled rather than dropped. True colour: a g, r, i → CIE XYZ matrix fitted to 129 galaxy spectra (leave-one-out error Δu′v′ 0.003, one just-noticeable difference) makes every galaxy a shade of white, 4,570 K (M104, M82) to 6,330 K (M101), matching their spectra to ≤ 0.007. Too faint for colour anyway: the brightest core (M104) reaches 0.08 cd/m², and at most 1.5% of any galaxy's area passes the 0.005 cd/m² cone threshold; no telescope can raise surface brightness. M51 on one grid from seven telescopes (FUV to 21 cm, ≤ 0.24″ from Gaia once stray SIP keywords are removed); Hubble's mosaic registered (0.84″ offset) and cross-calibrated against SDSS (−2/+9/+4% of the spectrum's prediction). Legacy Surveys DR9-north keeps only 36% of M51's z light. Plus an offline make-your-own-composite explorer. |
| G2 Meteorite map | **done** — [14-meteorite-landings/](14-meteorite-landings/). All 45,716 meteorites of the Meteoritical Society catalogue (NASA's 2013 snapshot; the old Socrata link is dead, the file moved to `data.nasa.gov/docs/legacy/`), on our own equal-area Equal Earth maps. Repaired first: 13,748 without usable coordinates (incl. 6,214 at a dummy 0°,0° and 223 Yamato meteorites whose zeroed latitude would put Antarctica in Kenya), one found on Mars, and every year before 1583 one too small (a Gregorian 1 January read back in the Julian calendar — all 8 match the arithmetic). 72% of located finds come from Antarctica (8% of the land); Oman yields 9,900 per million km², 150× the rest of the land. Against camera-network flux (58 falls ≥100 g per million km² per year), only 1 fall in ~1,500 is recorded, and recovery ∝ population density^0.59 (95%: 0.44–0.75) — 1 in 170 above 300 people/km², ~1 in 4,000 below 10; yet falls held at ~6/yr since 1900 while population grew. Evatt et al.'s iron table reproduced from the independent snapshot (falls 60/1,107 vs 60/1,097; Antarctica 0.74% vs 0.68%), and extended: individual finds are 40% metal vs 5% of falls; desert meteorites bought through dealers are 18% achondrites vs 4% where the find site was recorded. Plus an offline interactive map of all 31,893 located meteorites with filters. |
| H1 Sgr A* mass | **done** — [17-galactic-centre/](17-galactic-centre/). S2's orbit from both teams' raw measurements: 45 Keck and 145 NTT/VLT positions (1992–2018) plus 116 velocities — UCLA's Do et al. 2019 release (its page is gone; fetched from the Internet Archive) and Gillessen et al. 2017 via VizieR — fitted by MCMC with each team's reference frame free (one shared frame: χ² worse by 1,537). Each team's own analysis reproduced first: all 16 of Do et al.'s Newtonian posteriors within 0.15σ, their GR mass and distance within 0.04σ; Gillessen et al.'s VLT fit within 0.3σ on every number. Both teams together: M = 3.96 ± 0.04 ± 0.02 million M☉, R₀ = 7.93 ± 0.04 ± 0.03 kpc (Kepler by hand: 126.9 mas = 1,006 AU, P = 16.045 yr → 3.96 million suns). At GRAVITY's distance our banana gives 4.289 ± 0.010 million M☉ (GRAVITY: 4.297); their 4% larger R₀ is the known Keck–VLTI tension. Einstein's redshift detected: f = 0.82 ± 0.16 on Keck positions (Do et al.: 0.88), 0.70 ± 0.15 adding the VLT's, Newton excluded at 5σ; a Newtonian fit overstates M and R₀ by 5% and fakes Gillessen et al.'s 28 km/s "recession" (1.5 GM/ac = 17.5 km/s). A free precession fit claims 5× Einstein at 5σ — a false signal from 1990s speckle positions; adaptive-optics data alone cannot measure it. Pericentre 115 AU (16 light-hours, 1,470 horizon radii) at 7,590 km/s: > 5.5 × 10¹⁵ M☉/pc³, where the longest-lived dark cluster (Maoz 1998, reproduced) lasts 8 × 10⁵ yr. For I5: predicted shadow 51.2 ± 0.3 µas (EHT ring 51.8 ± 2.3). Next pericentre 5 June 2034. Plus an offline orbit explorer. |
| H2 21-cm rotation curve | **done** — [18-milky-way-rotation/](18-milky-way-rotation/). The HI4PI 21-cm survey: 25 cubes, of which only the \|v\| ≤ 350 km/s channels were fetched by HTTP byte range (3.8 of 6.6 GB). Terminal velocities from McClure-Griffiths & Dickey's two-error-function edge fit match the VLA survey at all 748 first-quadrant longitudes (median −0.24 km/s, rms 0.57), although HI4PI has a 16′ beam against the VLA's 1′. Their linear fits reproduced (v_circ 172.7 / 178.1 vs 172.5 / 177.7), as is the north–south asymmetry (−5.13 vs −5.2 km/s at \|l\| = 55–67°); an outward drift of the LSR explains none of it. Rotation curve from the tangent points with the LSRK convention handled explicitly (the Sun's angular speed, 30.32 km/s/kpc, from Sgr A*): 207–217 km/s at 2.6–3.6 kpc rising to 236–245 at 6–7.6 kpc, with ~10 km/s arm bumps; 63 of 1,442 longitudes with ambiguous edges dropped. The 3D orbits of 144 masers, computed from Reid et al.'s astrometry, lag rotation by 4.7 km/s and fall inward at 4.7 (published 4.3 ± 5.6, 6.0 ± 1.4); Perseus kinematic distances are 1.17× the parallaxes (up to 4.7×). Mass 1.06 × 10¹¹ M☉ inside the Sun's orbit; the Milky Way rotates faster than 89% of SPARC galaxies. Outer disc mapped face-on to 30 kpc after re-deriving LBH06's elliptical-orbit correction (with their constants α = 7.9, β = −1.06 vs published 8.67, −1.08; anticentre mismatch χ² 1137 → 0.6). HI scale length 4.05 kpc with the IAU constants (Kalberla & Dedes: 4.0 for co-rotating gas), 3.2 kpc with ours. Warp tilt W1 = 1.62 kpc at 20.5 kpc (LBH06: 1.74), maximum at φ ≈ 80–90°. The disc flares from 0.6 to 3.3 kpc thick. The plane cubes alone would have cut off the northern warp's top in 20–44% of cells beyond 16 kpc, so 7 cubes above the plane were added. |
| H3 Galactic archaeology | **done** — [19-galactic-archaeology/](19-galactic-archaeology/). Gaia DR3 only: 7.5 million rows from ESA's archive as small sky-patch queries, each all-sky sample checked against a server count. **GD-1** found blind among 3.8 million stars, with the ancient cluster M92 (same age, metallicity and distance) as the colour–magnitude filter and the stream's motion found in the data (up to 89σ): 3,464 members over 90° of sky; distance 7.8 → 11.6 kpc, gradient 0.0098 ± 0.0012 mag/deg (Koposov 2010: ~0.012), confirmed by the members' mean parallax (−0.004 ± 0.008 mas); gap G-20 78% deep (5.6σ), a 65% gap at −6.5°, the densest point at −14.5° (PWB2018's progenitor candidate: −13.5°); the spur 33 stars vs 2 in the mirror box; 65 of 66 spectroscopic members (Huang+2019) recovered. **Gaia-Sausage-Enceladus**: 74,938 halo stars within 2.5 kpc (0.41%); with no velocity cut at all, using Gaia's own RVS metallicities of cool giants, the halo's anisotropy jumps from β = 0.26–0.34 below [M/H] ≈ −2.2 to 0.74–0.83 above −2.0 (Belokurov+2018: 0.2–0.4 → ~0.9); the usual velocity cut finds only ~half of the round metal-poor halo. GSE median [M/H] −1.16 (Feuillet+2021: −1.17), 0.037 ± 0.010 dex poorer in α than the thick disc; its stars form the halo's blue sequence. Found on the way: half of GSP-Spec's "[M/H] < −1.66" stars are hot disc stars with wrong metallicities. **Phase spiral**: 2.5 million stars (Antoja+2018: 0.94 million); the snail shell in counts, v_R and v_φ once each star is compared with a phase-mixed disc; the disturbance dated to 250–560 Myr (median 360) in the McMillan 2017 potential, whose local disc matches the measured surface density (76 vs 68 ± 4 M☉/pc²); Antoja: 300–900 Myr. Plus an offline halo explorer. |
| H4 Sun's acceleration from quasars | **done** — [59-quasar-acceleration/](59-quasar-acceleration/). Gaia DR3's 1,215,942 five-parameter Gaia-CRF3 quasars (all three published counts reproduced from the archive), fitted with our own vector-spherical-harmonic code. The Sun accelerates at **5.05 ± 0.36 µas/yr = (2.32 ± 0.16) × 10⁻¹⁰ m/s²** towards (l, b) = (−1.1°, −3.3°), 3.4° from Sgr A*. Klioner et al. 2021's Table 2 is reproduced to 0.02σ in every component, error and correlation, as are their Gaussian widths (1.052/1.055/1.063), unit weight error and 3,783 outliers. Galactic constants predict 4.54–4.87 µas/yr (H1's R₀ + Sgr A*: 4.64 ± 0.09), all within 1.2σ. Implied: V_c = 238 ± 8 km/s, Galactic year 205 ± 7 Myr, 1.05 × 10¹¹ M☉ inside the Sun's orbit (H2: 1.06), and R₀ = 8.5 ± 0.6 kpc from the drift + Sgr A*'s motion alone. Gaia's own large-scale patterns (rms 6.0 µas/yr) exceed the signal (4.4); dipole-only fits point 12° east, as does Darling 2025's. Frame check: the 428,034 frame-rotator quasars show zero spin, but the spin runs with magnitude (−5 to +29 µas/yr). The Gaia-CRF3 spin check is reproduced exactly (1,607,289 sources) except the printed sign of ω_x. The drift depends on redshift at 4.0σ (Quaia: 3.50 ± 0.51 below z = 1.25, 6.27 ± 0.47 above; sharper with SDSS spectra), confirming Tsigkas Kouvelis et al. 2026's 2–3σ claim. Sky coverage (redshift-shuffle null), brightness and host galaxies are ruled out; the cause is open. Zero-parallax conditioning needs the zero point (naive: 8.2 µas/yr). Stars passing quasar cuts glide 140× faster; contamination < 0.5%. Plus an offline explorer. Spec corrections above: V₀ source, ~0.6 GB, watch-outs. |
| H5 Multiwavelength Milky Way | not started |
| H6 3D dust fly-through | **done** — [60-dust-flythrough/](60-dust-flythrough/). The Edenhofer et al. (2024) map (main + 2 kpc maps, 7.4 GB) sampled exactly as the authors' tools do (identical to `interp2box.py`, 6e-8 from `dustmaps`). Keeper: a WebGL 2 viewer (no libraries; three nested volumes, 2.5-10 pc cells; Glow / Starlight / Solid modes, each stated) and a 100-s movie recorded from it on the GPU. Their Fig. 4, read back through its colour bar, matches ours (r = 0.995) but is 4-8% lower, rising with distance: 12 posterior samples fetched by byte range show that the figure averages log-interpolated samples while the released mean map (and `interp2box.py`, `dustmaps`) interpolates the mean, which Hölder's inequality makes larger by exactly that much. Local Bubble: O'Neill et al. (2024)'s recipe re-run on 786,432 sightlines, 97.5% within 1 pc of their model (median 169 pc, thickness 35 pc, 6.0e5 M☉, interior = a 166-pc sphere); Zucker et al. (2022)'s expansion refit from their tracebacks (14.5 Myr, 166 ± 6 pc, 6.7 km/s); young clusters within 250 pc lie within 30 pc of the wall 51% of the time (random 35 ± 8%, 60-300 Myr 29%). Radcliffe Wave: unrolled in the dust (no fit); Konietzka et al.'s model fitted to Zucker et al.'s 326 sightlines with map distances: 2.95 kpc long, up to 203 pc, mean wavelength 2.16 kpc (2.73 kpc, 169 pc with the published distances); Alves et al.'s printed eq. (1) + Table 2 (δ = 6.68) never exceed 61 pc, and a refit with their priors gives δ = 0.30; broad priors prefer other shapes (the smooth wave is partly the priors). Clouds 2-4% farther than Zucker et al. within 700 pc. Vergely et al. (2022): integrated extinction 0.998, but radial elongation 2.2-3.0 vs 0.84-1.27. Dereddening tool validated on 1,312 clusters (0.986 once R_V 2.76 vs 3.1 is accounted for) and SFD (slope 1.045, +0.008 floor). Found: `dustmaps`' integrated mode is half a voxel off, `interp2lbd.py` skips radial interpolation, posterior voxels correlate +0.3-0.5. |
| I1 M–σ relation | **done** — [32-m-sigma-relation/](32-m-sigma-relation/). Three compilations rebuilt from their own tables (van den Bosch 2016 FITS; Kormendy & Ho 2013 Tables 2-3 from the arXiv source, Table 3's fit flags read off the colours of its PostScript; McConnell & Ma 2013 Table A1) and refitted: every published slope back within 0.11σ (4.36 vs 4.38, 5.66 vs 5.64, 5.34 vs 5.35; van den Bosch needs his true-σ modelling, a plain likelihood dilutes the slope to 5.20, shown by simulation). The slopes differ by sample: from 5.66 to 4.36, dropping pseudobulges gives −0.63 and four new extreme-σ galaxies −0.40. Ellipticals + classical bulges: slope 4.37 ± 0.31, scatter 0.29 dex; pseudobulges sit 0.60 dex lower with no trend of their own. Sgr A* (H1 mass) 0.68 dex below the classical relation, M87* (EHT) 0.40 above. Black hole = 0.49% of the bulge (Kormendy & Ho method; dynamical bulge masses give 0.28%). King's feedback M_σ predicts 3.0×10⁸ M☉ at 200 km/s with no free parameter (measured 3.15×10⁸). ATLAS³ᴰ test of selection: black hole hosts have σ 0.050 ± 0.014 dex higher than other galaxies of their luminosity (p = 0.0006); a labelled simulation shows unresolved spheres of influence flatten detection-only slopes. Reverberation ⟨f⟩ = 2.2-9.8 depending on the quiet galaxies matched. Interactive refit-it-yourself page. |
| I2 Reverberation mapping | **done** — [33-reverberation-mapping/](33-reverberation-mapping/). Our own interpolated cross-correlation code (PyCCF's rules; equal to a literal loop to 10⁻¹⁶) on three data sets. NGC 5548 in 2014 (AGN STORM, 133 MDM nights): Hβ lag 4.16 (−0.36/+0.37) d (Pei et al. 2017: 4.17 ± 0.36), so the broad-line gas sits 721 AU from the black hole; virial product 1.49 ± 0.48 × 10⁷ M☉ (Pei: 1.49 ± 0.49), mass 6.3 × 10⁷ M☉ with our f. The "BLR holiday" measured against HST's far-UV continuum: flux lost C IV 17.7%, Si IV 23.7%, He II 22.9%, Lyα 10.0%, Hβ 6.2% (Goad et al. 2016, Pei et al.: 18/23/21/9/6%). The disk echoes too: 12 continuum bands lag the far-UV as λ^4/3 with α = 0.427 ± 0.015 d (Fausnaugh et al. 2016: 0.42 ± 0.02), a disk 2.8× larger than thin-disk theory whose optical light comes from as far out as the inner broad-line gas. AGN Watch 1988–2008: 15 of 16 published seasonal lags reproduced within 1σ (median 0.08 d; only without detrending); NGC 5548's own R ∝ L^(0.52 ± 0.20), the photoionization 0.5; Bentz et al. 2007's 0.73 comes from counting the common host-starlight error per season (their published numbers through our code: 0.71); in 2014 the lag was 5.0× shorter than its own trend. Virial test: 29 line/season pairs over 25 years give one mass to 0.21 dex (slope −0.38 ± 0.07 vs gravity's −0.5). R–L: LINMIX returns Bentz et al. 2013 exactly (1.527, 0.533, 0.19); the database's 84 host-corrected lags of 47 AGN give slope 0.528 ± 0.034, scatter 0.21 dex; SDSS-RM Hβ, Mg II and C IV reproduce Shen et al. 2024 to 0.01 (slopes 0.42, 0.31, 0.32; C IV scatter 0.51 dex), SDSS-RM Hβ 0.08 dex below the local relation. 444 SDSS-RM lags remeasured from raw light curves with alias weighting: median ratio to the published ICCF lags 0.995–1.010, 95–96% within 1σ; false-positive rate ~7% (118 positive vs 8 negative lags at r > 0.4). Virial factor: ⟨log f⟩ = 0.626 ± 0.075 from 30 dynamical models, scatter 0.33 dex (Shen: 0.62, 0.31); 0.75 ± 0.12 from M–σ with van den Bosch's 160 galaxies (Woo et al. 2015 with Kormendy & Ho's: 0.65). NGC 5548 sits 0.39 dex below M–σ. Found on the way: the LINMIX Python port ignores its seed when run serially (fixed). Plus an offline find-the-echo page. |
| I3 Quasars | **done** — [34-quasars/](34-quasars/). SDSS DR16Q (750,414 quasars) with Wu & Shen 2022's measurements (1.39 GB, read column by column in one streaming pass), the DR7 uniform sample (57,959) and 114 quasars beyond z = 5.7 (Farina et al. 2022). Richards et al. 2006's luminosity function rebuilt from their own 15,343 quasars (93% of 101 bins within 0.05 dex; one published bin contradicts its own stated limit); our DR7 one matches Shen & Kelly 2012 to −0.0003 ± 0.009 dex once their omitted 5% image factor is noted. Luminous quasars (M_i < −27.6) peak at z = 2.59 ± 0.01 (2.50–2.73 across analysis choices), 2.5 Gyr after the Big Bang, 60 per Gpc³, a billion years before cosmic star formation. Our own composite of 4,446 spectra matches Vanden Berk et al. 2001 to 1–2% over 1,250–3,500 Å (luminous quasars: α_ν = −0.48 vs −0.44) and shows the Baldwin effect; its colours, moved through the SDSS filters, enter the stellar locus at z = 2.70, exactly where the selection function bottoms out (5%), and its K-correction matches Richards et al.'s to 0.04 mag. 707,947 masses (recomputed identically to the catalogue): Eddington ratios lognormal, median 0.13, 2.8% above the limit. C IV masses grow +0.17 dex per 1000 km/s of outflow blueshift at fixed luminosity; Coatman et al. 2017's correction over-corrects these fainter quasars, our cross-validated recalibration removes the trend and an artificial 0.24 dex mass drop at z = 2. J0313−1806 (671 Myr, 1.6 × 10⁹ M☉) needs a seed ≥ 1.7 × 10⁴ M☉ even growing at the Eddington limit from z = 30 (Wang et al. 2021: 10⁴–10⁵); 84% of the z > 5.7 quasars need more than 100 M☉. Plus an offline quasar time machine page. Download ~2.9 GB (spec: ~1–2 GB). |
| I4 X-ray binaries | **done** — [35-x-ray-binaries/](35-x-ray-binaries/). All 20 dynamically confirmed black holes of the online BlackCAT (Sept 2025; 73 transients), parsed with every value's paper, plus Cyg X-1, LMC X-1/X-3, M33 X-7 and Gaia BH1–3. Mass functions recomputed from each listed P and K₂: 10 exceed 3.2 M☉ (13 exceed 2.9), V404 Cyg 6.068 ± 0.044 (Casares & Charles: 6.08 ± 0.06); Cyg X-1 from its orbit 21.2 ± 2.3 M☉ (Miller-Jones et al.: 21.2 ± 2.2). BlackCAT inconsistencies found: GRO J1655-40 pairs a 2008 K₂ with a 1999 f(M) (4.3σ; mass 7.8 vs listed 6.0), V4641 Sgr's q contradicts its cited paper (override → 6.8 M☉ vs MacDonald's 6.4), GS 2000+251's mass range does not follow from its listed i. Farr et al. 2011's hierarchical mass-gap test replicated with nested sampling (exponential and Gaussian M₁% quantiles to 0.1 M☉; power-law 10% quantile 3.9 vs 4.3); Kreidberg et al. 2012's inputs remove the gap (M₁% 1.5/2.4/3.5); today's BlackCAT keeps it (M₁% > 5.5 power law, > 3.7 Gaussian at 90%), while CatNoir 2026 inputs for the same 19 systems give > 2.54 (Rocha et al. 2026: 2.56). 12 of 19 black holes cannot reach 5 M☉ at any inclination. Against 282 GWOSC mergers: median 28.6 M☉ as detected vs 8.3 for X-ray binaries, 10.6 after 1/V_max weighting (heavier mergers are louder). MAXI/BAT: MAXI J1820+070 hard→soft 2018-07-07 (Homan et al.: July 3–9), peak 9.9% L_Edd (2–50 keV), soft→hard at 1.3%; GX 339-4 11 outbursts since 2009, 6 full loops, hysteresis 4–18×. Gaia BH1 from DR3 alone: a₀ 3.01 mas, i 121.2°, M₂ 13.2 +2.9 −2.2 (El-Badry et al.: 12.8 ± 2.0); BH2 f_ast 6.73 (6.72). Spec corrected: BH1's 9.6 M☉ needs follow-up velocities. Offline explorer page. |
| I5 M87* shadow | **done** — [36-m87-shadow/](36-m87-shadow/). The EHT's public calibrated visibilities (M87 2017 and 2018, Sgr A* 2017; 11 MB) analysed from scratch: our reader and scan averaging reproduce Paper VI's visibility counts exactly on all eight 2017 data sets, and our closure phases equal ehtim's to 10⁻⁶°. The amplitudes alone show a ring: first null at 3.4 Gλ → 46 µas. Twelve models ranked by Bayesian evidence on 11 April: lopsided rings with faint extra emission beat disks and blobs by thousands of log-units. M87's ring from 52 nested-sampling fits, each checked for missed modes (6 caught and fixed): 43.8 µas (night-to-night rms 1.1) for a blurred thin ring; Paper VI's crescent model 43.8 ± 0.2, within 0.6 µas of Paper VI on 7 of 8 data sets (April 11 at 227 GHz: 45.0 vs 41.6, where Paper VI's own two variants disagree); amplitude nulls alone 43.9. Between April 5–6 and 10–11 the ring grew 1–2 µas and its bright side swung from ≈140–160° to 180–200° (Paper VI: the same). With the EHT's GRMHD calibration (α = 11.56): θ_g = 3.77 (+0.49/−0.39) µas, M = 6.4 (+0.9/−0.7) × 10⁹ M☉ (EHT: 3.8 ± 0.4, 6.5 ± 0.7); 7.1 ± 0.4 if the ring were the shadow's edge. The ring sides with stellar dynamics (Gebhardt et al. 2011: 0.2σ) against gas dynamics (Walsh et al. 2013: 2.6σ; its mass would make a 24 µas ring). 2018: 44.2 µas, bright side moved to 202–216° (EHT: 209–222°). Images from our own RML imager and from the EHT's 2019 eht-imaging script re-run here (its image fits our closures at χ²/N = 1.00) both show the ring (39–42 µas peak to peak); simulated disks, blobs and doubles do not come back as rings. Sgr A*: closure-only fits are degenerate (a ring or a tiny crescent in a halo: 56 vs 36 µas), so the ring is read from the amplitude nulls, with scattering and Paper IV's noise budget: 53.9 ± 1.3 µas (EHT: 51.8 ± 2.3) = 10.95 ± 0.27 θ_g of H1's S2-orbit mass, 1.4 ± 2.5% above the ring the EHT calibration predicts for it (53.1 µas). ehtim's scattering class fails on numpy 2.5, so it was reimplemented (agrees to 6 × 10⁻⁷). Plus an offline fit-the-ring-by-hand explorer. |
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
| I18 Giant radio galaxies | **done** — [21-giant-radio-galaxies/](21-giant-radio-galaxies/). All 11,585 known giants (Mostert et al. 2024) plus LoTSS DR2 images of 17 we measure ourselves. Porphyrion: 13.6′ at PA 26° (published 27°) → 6.54 Mpc projected (Oei et al. 2024: 6.43); for random orientations, most likely 6.8 Mpc in 3D (published 6.8). Oei et al.'s length-distribution model rebuilt and refitted on their 1,473 giants (a 61-million-point grid): tail index ξ = −3.45 ± 0.57 vs their −3.5 ± 0.5 (−2.72 if selection is ignored); 12 simulated surveys confirm the method (mode unbiased). The record claims: four of the five longest catalogue giants (7.98–5.66 Mpc) rest on detached sources; the host's own emission spans 2.2–3.5 Mpc. Catalogue lengths reproduced to 0.014%; with D1's H0, 19% would no longer be giants; 2.9% of SDSS-observed hosts carry redshifts >10% off. 3C 236's lobes: spectra steepen toward the host (JP ages ~30 → 80 Myr); our CI model on Shulevski et al.'s fluxes gives 136/114 Myr vs their 129/117. Black holes (M–σ, 2,053 hosts): median 0.91 × 10⁹ M☉, the same as ordinary radio galaxies at equal z and radio power (+0.016 ± 0.012 dex), and unrelated to length. Environments: Oei et al.'s BORG densities recomputed (giants 4.73 vs 3.84, p = 3 × 10⁻⁶), but equal to radio galaxies above 10²⁴ W/Hz (4.44, p = 0.60); our SDSS neighbour counts agree (0.95 ± 0.07 vs same-z radio galaxies; 1.35 ± 0.10 vs radio-quiet look-alikes). LoTSS DR3 (2026, 88% of the northern sky) adds ten giants that DR2's public images miss: NGC 2300's, the largest on the sky (a four-cutout mosaic near the pole: 117′, published 2.2°), NGC 6251 (1.90 Mpc), DA 240, J1129+2733 (3.66 Mpc) and others. J1526+2620, claimed at 4.56 Mpc, shows only a 2.44-Mpc arm. The 17 DR2 giants re-measured on DR3 agree in length to 0.1% and in flux to 0.2% (median). Plus an offline explorer of all 20 giants with a measuring tool. |
| J1 Velocity fields | **done** — [25-velocity-fields/](25-velocity-fields/). THINGS 21-cm cubes of NGC 3198, NGC 2403, NGC 2903 and DDO 154 (2.1 GB from MPIA): 3.1 million line profiles fitted with our own vectorized Gauss–Hermite/Gaussian solver, then tilted rings (our ROTCUR, de Blok et al.'s sequence) from our own starting guesses. NGC 3198's cube labels observed frequencies: read naively it is 24.5 km/s off; via AIPS's ALTRVAL every channel matches the THINGS global spectra to 0.1 km/s. de Blok et al. 2008's curves exist only as figures, so their points were read back from the SuperMongo vector PostScript (ring spacing recovered exactly): ours agree to 0.95 / 2.60 / 3.15 / 0.92 km/s rms over 512 rings, 99–100% within 1σ; V_sys within 0.0–2.3 km/s, centres 0.3–9″. Five velocity recipes compared: only the Gauss–Hermite *peak* reproduces their NGC 2403 curve (2.6 vs 5.3 km/s rms for the centre parameter), so that is what they used; THINGS' first-moment maps drag inner rises down by 10–12 km/s. NGC 2903's bar needs the geometry frozen inside 110″ (as they did; −5 km/s otherwise). Harmonic decomposition (RESWRI-style): median non-circular A_r 4.74 / 4.46 / 8.98 / 2.04 km/s vs Trachternach et al.'s 4.49 / 4.03 / 6.10 / 1.61; the bar's inner kpc 14.4 vs 13.6. HI fluxes within 0.1–3.8% of Walter et al. DDO 154: cored halo R_c 1.25 ± 0.09 kpc, ρ₀ 36 ± 4 ×10⁻³ M☉/pc³ (Oh et al. 2011: 1.30, 34.5); NFW only with c = 5.8 ± 0.6 vs ΛCDM's 11; inner slope α = −0.40 ± 0.13 (Oh's −0.29 ± 0.15 read from their figure; derivative choices move α by ~0.3). Labelled simulation: a ΛCDM cusp built into a model cube, smoothed with the real beam and run through the identical pipeline reads α = −0.77 and still spins 1.55× faster than DDO 154 inside 1 kpc (χ² = 41 for 7 points): beam smearing cannot hide a cusp. MaNGA DR17: 10,461 galaxies streamed (only the 38% of each MAPS file needed, 8 GB in an hour), kinematic axes of stars and H-alpha by pafit's method (identical results, 8× faster) plus the sense of rotation: among 5,226 regular rotators 7.4% have misaligned gas (Zhou et al. 2022: 6.6%), 2.3% counter-rotating (120); blue 2.1% (Zhou: 2.4%), red 29% (ATLAS3D: 36 ± 5%); without vetting 22.7%. Asymmetric drift in three discs: stars lag gas by 7–18 km/s, 4–7× less than the Jeans equation predicts from MaNGA's dispersions, which sit at its ~70 km/s resolution. Plus an offline tilt-a-disc explorer with spectra. |
| J2 Galaxy morphology | **done** — [23-galaxy-morphology/](23-galaxy-morphology/). Galaxy Zoo 2's votes for 239,695 galaxies (Hart et al. 2016 debiasing) joined to SDSS DR17, MPA-JHU stellar masses and D2's density measure (its 134,975 tracers and red/blue divider reproduced exactly). Distant galaxies look smoother: for one kind of galaxy the raw 'features' vote falls from 0.69 to 0.28 between z = 0.026 and 0.086; in a volume-limited flatness test Hart's correction still drifts 4–9 points (Willett's up to 11, raw votes 27). Published results rebuilt with their own cuts: red face-on disky spirals 5.9 ± 0.4% (Masters+10: 6%), blue ellipticals 7.0 ± 0.4% (Schawinski+09: 5.7 ± 0.4%), bar fraction 28.7% (Masters+11: 29.4%, once the vote threshold matches their half-complete data), two-armed spirals 62.0 ± 0.4% (Hart+16: 62.1 ± 0.4%) with many-armed ones bluer at fixed mass, and 11 of 12 numbers of Willett+13's comparison with an expert (e.g. strong/intermediate/weak bars with p_bar > 0.8: 78/40/9%, exact). Bars: 39% of bright face-on disks, most common at intermediate colour; at 10^10 M☉ red disks are barred 2.5× as often as blue. Environment at fixed mass: the red fraction rises 6–26 points from the sparsest to the densest fifth, 75–90% of it within each shape (an exact decomposition, robust to dropping edge-on disks), Bamford+09's result with D2's density. ML extension: a 0.39-million-weight CNN trained on 9,000 volunteer images (pulled one by one out of the 3.35 GB Zenodo zip by byte range) matches 99.1% of clear smooth-or-features majorities on 2,063 unseen galaxies (Zoobot on deeper DESI images: 95.6%, seeing +0.11 more disks), and inherits the volunteers' redshift bias, as Zoobot does. Plus a tuning-fork poster chosen purely by vote fractions and an offline 'classify like a volunteer' explorer. |
| J3 Galaxy life cycle | **done** — [24-galaxy-life-cycle/](24-galaxy-life-cycle/). 656,838 SDSS DR17 Main Galaxy Sample galaxies with MPA-JHU masses and SFRs (MPA-JHU tables fetched by specObjID range and joined locally), cross-checked against GSWLC-2; completeness and area from the survey's own target list (776,234 targets, 92.9% observed) and tiling polygons (8,031 deg², DR7 published 8,032). Main-sequence ridge line (Renzini & Peng's definition, 1/V_max-weighted): slope 0.69 ± 0.02 over 10⁹–10¹⁰·⁵ M☉ (GSWLC 0.75; R&P 0.76), 0.62 over 10⁹–10¹¹ because the ridge bends where the star-forming peak fades. Green valley 0.86–0.99 dex below the sequence at every mass from 10⁹·² to 10¹¹; half of galaxies quenched above 10¹⁰·⁴⁰ M☉ (GSWLC has no valley at all; same cut gives 10¹¹·⁰⁵). Peng et al.'s mass/environment separability confirmed: environment threshold log p₁ = 1.83 ± 0.03 (Peng: 1.84), mass dependence of the environmental effect −0.07 ± 0.05. "AGN crowd the green valley" holds with MPA-JHU SFRs (33% vs 19% on the sequence) but vanishes with GSWLC's UV-to-IR SFRs (19% vs 19%). Local SFR density 0.0145 ± 0.0002 M☉/yr/Mpc³ (Kroupa, likelihood means; medians 30% lower); Salpeter-converted estimates 0.015–0.022 bracket Madau & Dickinson's 0.0193; stellar mass density log 8.588 (Li & White: 8.59). 67 objects with log sSFR > −8 (one claiming 7,500 M☉/yr) dropped. Madau & Dickinson's Table 1 parsed from their LaTeX and refitted: peak z = 1.72 (1.63–1.81), 3.7 Gyr after the Big Bang; half of today's stars by z = 1.24; their own curve, through our code, reproduces every milestone they quote (z 1.86, 3.45 Gyr, 5.76 × 10⁸ M☉/Mpc³, 0.81% before reionization) except the mean stellar age (8.04 vs 8.3 Gyr), yet fits their table worse (χ² 64 vs our 40). Direct stellar-mass counts sit 0.11–0.18 dex below the integrated history at all epochs. Plus a rule-picked SDSS image gallery, a quenched-galaxies cosmic-web slice and an offline explorer. |
| J4 Laniakea flows | **done** — [22-laniakea-flows/](22-laniakea-flows/). Cosmicflows-4 (55,877 galaxy distances) against the 2MRS mass map (44,599 galaxies). Two problems in CF4 as published (VizieR/journal tables) found and documented: its Watkins–Feldman velocity columns use log10 instead of ln (2.30× too small; reproduced to 0.25 km/s; their "300 km/s" scatter is below the distance errors alone, 508), and its 34,045 SDSS Fundamental Plane moduli carry the wrong sign (DM = DM_z + 4.98 η, r = 0.9975, while geometry fixes η = log10(d_z/d_H): Howlett's sizes are reproduced at the group distance to 0.0002 dex). Corrected and re-levelled (−0.055 ± 0.016 mag), our fits reproduce Nusser & Tully's 2026 bulk-flow analysis to 0.33σ (1.19σ as published). H0 = 74.37 ± 0.08 ± 0.75 (zero point) (paper 74.6; table as published 74.99); Planck's H0 would need a 0.22 mag shorter ladder. 2MRS flux dipole 11.6° from our CMB motion at 60 h⁻¹ Mpc and 20.1° at 130 (Erdoğdu et al. 2006: 12 ± 7°, 21 ± 8°), β_L = 0.42 ± 0.04 (0.40 ± 0.09). Real-space 2MRS reconstruction vs 14,837 CF4 groups: β* = 0.407 ± 0.009 (Carrick et al. 2015: 0.431 ± 0.021), external flow 135 km/s toward (310°, −24°); predicted Local Group motion 660 km/s, 1.3° from the CMB's 620; measured vs predicted flow maps correlate r = 0.68. Bulk flows: 572 ± 83 km/s at 120–160 h⁻¹ Mpc (2.9σ; N&T 628 ± 82), all shells 2.4σ jointly (N&T 2.4σ), driven by TF and 6dFGS (SNe consistent); spheres reach ~300 km/s, inside ΛCDM's 68–97.7% band for their real window. Basins of attraction by watershed of the 2MRS potential: at 6 h⁻¹ Mpc smoothing our basin spans 12,600 km/s and 1.2 × 10¹⁷ M☉ around an attractor beside Centaurus (Laniakea: 12,000 km/s, 10¹⁷); it grows with smoothing and, with the external pull, drains to Shapley (cf. Valade et al. 2024). Plus an offline 3D explorer. |
| J5 Einstein rings | **done** — [26-einstein-rings/](26-einstein-rings/). All 85 SLACS definite lenses: SDSS DR17 spectra of the 131 candidates Hubble imaged, and Hubble cutouts (MAST mosaics) of 77. A blind search of the spectra recovers the background galaxy's redshift in 130 of 131 (after masking the lens's own emission lines); confirmation rises from 25% for predicted rings < 0.8″ to 86% above 1.2″. Ring masses (2.6 × 10¹⁰–5.4 × 10¹¹ M☉ inside 1.1–7.2 kpc) reproduce SLACS IX to 0.0006 dex. Our own lens modeller (SIE + Sérsic source + multi-Gaussian lens light; deflections match lenstronomy to 10⁻⁹″; an injected lens recovered exactly) re-measures 67 published Einstein radii: median −0.7%, robust scatter 2.9%, 2.0% for the 42 fits that pass their own checks (SLACS adopted 2%), but 16 off by >10%, mostly disc and spiral lenses; external shear moves radii by 2.5% lens to lens. Isothermal test on Bolton et al.'s 53 lenses: σ★ = 0.948 ± 0.008 σ_SIE in the fibre, 1.016 ± 0.008 at R_e/8 (published 0.948, 1.019), intrinsic scatter 7.5% (published 7.5%); 1.008 ± 0.010 with our own radii. Lensing + Jeans dynamics: ⟨γ′⟩ = 2.097 ± 0.029, intrinsic scatter 0.17 (Auger et al. 2010: 2.078 ± 0.027, 0.16); lens by lens 0.020 rms from their table. Inside the ring stars are 43% of the mass with a Chabrier IMF, 76% with Salpeter; stars + NFW halo need log α = +0.024 ± 0.016 relative to Salpeter (Treu et al. 2010: 0.00 ± 0.03), 17σ above Chabrier, rising with σ; the two components sum to slope 2.15 ± 0.20, the bulge–halo conspiracy. Plus an offline explorer for dragging a galaxy behind a real lens. |
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
| L1 Supernova nickel | **done** — [27-supernova-nickel/](27-supernova-nickel/). Bolometric light curves integrated from UV, optical and infrared photometry of SN 2011fe, SN 2014J and SN 2023ixf, plus SN 1987A's published UVOIR. Nickel-56: SN 2011fe 0.58 ± 0.03 Msun from Arnett's model, which also gives the gamma-ray escape time T0 = 33.7 ± 1.6 d (published 34.5); SN 2014J 0.60 ± 0.09 from its light and 0.65 ± 0.17 from INTEGRAL's cobalt-56 lines (Churazov et al. 2014: 0.6 ± 0.1); SN 2023ixf 0.064 ± 0.003 (T0 = 210 ± 5 d); SN 1987A 0.066 ± 0.002 (textbook ~0.07). Arnett's rule for 90 low-dust ZTF Type Ia: median 0.52 Msun, and more nickel goes with slower fading (r = −0.44). After day ~150 part of SN 2011fe's light goes missing (far infrared or escaping positrons). Plus an offline page for building a light curve from nickel. |
| L2 Supernova census | **done** — [30-supernova-census/](30-supernova-census/). The ZTF Bright Transient Survey catalogue (20,599 candidates; 5,168 in the statistical sample, June 2018 - October 2025, 95% classified), its packed quality bits decoded and its absolute magnitudes reproduced to 0.0007 mag; Cosmicflows-4 distances for the 154 nearest (SN 2023ixf: −15.7 → −17.3). What the survey sees, Ia/Ibc/II = 77/5/18%, matches Li et al. 2011's magnitude-limited prediction (79/4/17); per unit volume (1/V_max, M < −14) the mix is 30/16/54% (LOSS: 24/19/57). Perley et al.'s rates re-run on their window: Ia 2.20 (2.03-2.41) × 10⁴ Gpc⁻³ yr⁻¹ vs 2.35 ± 0.24, CC 9.4 vs 10.1, recovery 414/691 vs 412/690. Season by season the survey is stable from June 2019 (5% rms) but its first season found 24% fewer; 2019-25: Ia 3.14 × 10⁴ (20-30% above PTF/ASAS-SN: the sky coverage is borrowed from 2018-20 because IRSA's exposure log was down), CC 9.4 × 10⁴; the coverage-free ratio CC/Ia = 3.00 ± 0.3 gives 7.3-7.7 × 10⁴ with PTF's or ASAS-SN's Ia rate, i.e. 7.8-8.7 core-collapse SNe per 1,000 M☉ of stars formed (Kroupa, 8-25 M☉: 8.8). Hosts (SDSS + MPA-JHU, 1,050): per unit stellar mass Type Ia rise 16× from quenched to strongly star-forming hosts (Mannucci et al. 2005: ~20×), core collapse 300×; Perley et al.'s host table reproduced (red sequence 32.8 / 12.3% vs 31 / 11%); six of the eight 'dead' core-collapse hosts are spirals whose bulge filled SDSS's 3″ fibre. Type II plateaus from ALeRCE light curves: 65 complete; the fall begins at 87 ± 16 d (Anderson et al. 2014: 83.7 ± 16.7), ⁵⁶Ni 0.065 M☉ median, E ≈ 1 × 10⁵¹ erg and M_ej ≈ 7 M☉ for a 500 R☉ star (12 M☉ for 300 R☉: Goldberg et al.'s scalings leave the radius free). Red supergiant problem from the tables of Smartt (2015) and Davies & Beasor (2018), parsed from arXiv: Smartt's '13 expected, p = 3 × 10⁻⁵' reproduced (13.6, 1.8 × 10⁻⁵); his likelihood gives M_hi = 15-16 M☉ (published 16.5), Davies & Beasor's late-M luminosities 19.9 (< 26.4; published 19.0, < 27); luminosity limit log L = 5.27 (theirs 5.28), 5.5 disfavoured at 1.5σ; a 24-SN sample would show nothing brighter than seen 11% of the time (1.2σ), and Beasor et al.'s extra scatter does not help (7%). Plus an offline explorer: click any supernova for its light curve and host. |
| L3 Remnant expansion | **done** — [28-supernova-remnants/](28-supernova-remnants/). Palomar plates (POSS-I 1951, POSS-II 1991, from the DSS) and Pan-STARRS (~2012) registered on Gaia DR3 stars wound back to each plate's date (0.15″ / 0.10″ / 0.015″ rms; bright stars show the photographic magnitude equation), then 810 Crab filaments tracked by masked correlation over 60 years: they converge in 1106.0 ± 7.8 CE (Martin et al. 2025: 1105.5 ± 0.5; Nugent 1998: 1130 ± 16); the other red pairs give 1105 (1951→1991) and 1098 (1991→2012); closest approach 1118 (theirs 1140.5). Today's speed is 1.059× the mean since 1054 (their v/v_h = 1.06): the pulsar wind nebula has pushed the filaments; median 1,330 km/s at 2 kpc. The Crab pulsar's Gaia DR3 motion (111 km/s) runs back to 2.1″ from our expansion centre in 1054. Cas A in 15 Chandra epochs (2000–2019; the archive's 0.7 MB 'centre images' instead of 520 MB event lists), Poisson stretch fits about Thorstensen's centre: forward shock 0.185%/yr, m = 0.62, ~5,000 km/s (Vink et al. 2022: 0.218, m = 0.73), 13 of 18 sectors within 0.03%/yr of theirs (the west is contaminated by debris in all-energy images); inner ring follows their reverse-shock pattern (r = 0.89). Neutron star 7.0″ from the centre at PA 169° → 340 km/s (Fesen et al. 2006: 7″, 169°, ~350); its direct proper motion is not measurable here (pointing/asymmetric-expansion degeneracy, shown). Cas A's optical shell on survey plates: a documented failure (every pair peaks at zero stretch; the same code recovers the Crab's 1.070 vs 1.071 predicted). Green's 310 remnants on the sky and nine historical ones against free expansion and Sedov–Taylor. ATNF v2.8.1: young pulsars 353 ± 43 km/s mean 2-D speed (Hobbs et al. 2005: 307 ± 47); 8 of 9 secure remnant associations move within 35° of straight outward (p = 1.5 × 10⁻⁵). Plus a blink comparator with a 'rewind' slider and a 15-epoch Chandra movie. |
| L4 Origin of the elements | **done** — [29-element-origins/](29-element-origins/). GALAH DR4: 530,403 clean field stars (of 917,588; 178,243 giants, 277,323 dwarfs) with their orbits, plus APOGEE DR17 (VizieR) and H3's Gaia-Sausage-Enceladus boxes rerun in H3's own code (77,205 / 12,343, exact). α-knee from both disc sequences fitted together: thick disc at [Fe/H] = −0.84 (GALAH giants; dwarfs −0.89, ridge-line fits −1.12/−0.81, APOGEE giants −0.49: the disc's knee is survey- and method-dependent); GSE at −1.33 (−1.43…−1.21; Koppelman box −1.43; Mackereth+19: ≈−1.2), 0.5 dex lower within the same survey: GSE formed stars more slowly. Two-process model (Weinberg+19) for 26 elements: Griffith+22's recipe on DR4 dwarfs reproduces 13 of their 16 DR3 prompt shares within 0.07 (Fe 0.51 vs 0.51, Ba 0.30 vs 0.31); giants give Fe 0.50, Mn 0.23, Ba 0.24 and **europium 0.98 ± 0.01** — in the disc, Eu arrived as promptly as Mg (giants vs dwarfs differ by up to 0.19 for O, Ti, Y, Cu). GALAH DR4's own Eu is unreliable (r = 0.40 with DR3; prompt share 1.27): Kane et al. 2026's VAC used instead. GSE − disc at fixed [Fe/H]: [Eu/Mg] +0.24 ± 0.02, Al −0.33, Na −0.18, Mn and Cr identical. JINAbase: [Eu/Fe] spreads 0.6 dex below [Fe/H] −2.5 (Mg 0.2); [Ba/Eu] rises from the pure-r −0.80 to solar. AT2017gfo: Villar+17's 646 points through 65 SVO filters, their three-component model re-implemented (Arnett integral done exactly; a trapezoid overshoots 15% late on) and sampled by nested sampling: blue 0.020, purple 0.045, red 0.013, total 0.077 ± 0.002 M☉ (Villar: 0.020/0.047/0.011/0.078); a 1% second mode swaps purple and red. Europium budget: 16 ± 3 M☉ of r-process Eu in the Galaxy vs 9.3 × 10⁻⁵ M☉ per GW170817-like merger → supply/demand 0.14–4.7 over GWTC-4.0's rate range (0.82 at its middle); Villar+17's 2017 "factor-10 surplus" becomes a shortfall at today's 35× lower rate. Plus a periodic table of origins coloured by the measurements and an offline element explorer. |
| L5 Type Ia survivors | **done** — [31-hypervelocity-survivors/](31-hypervelocity-survivors/). Gaia DR3 via the ARI Heidelberg mirror (ESA overloaded, VizieR down; whole-sky cuts run patch by patch on `gaia_source_lite`). 13,538 stars move faster than 600 km/s across the sky at face value; for 76% a parallax just 1σ larger makes them ordinary halo stars, and only 5 have a parallax < −3σ that would pass (wild errors are rare). With the Lindegren zero point and Maíz Apellániz et al.'s realistic errors, 22 stay robust and 5 must be unbound whatever their radial velocity: D6-1, -2, -3, the bound eclipsing binary J1903−0023 (Bhat et al. 2025: parallax misleading) and a halo hypervelocity star (Carballo-Bello et al. 2025). El-Badry et al.'s ADQL query returns the same 25 sources today; Shen et al.'s four DR2 impostors show bad DR3 astrometry (parallaxes +3.3–5.7σ with RUWE 1.4–3.8; O3's 285 mas/yr was a blend 0.55″ from a brighter star). 3D velocities of all 17 proposed survivors reproduce El-Badry et al.'s Table 3 to ~1%, but only with uncorrected parallaxes (with the zero point, D6-1 launches at 2,152 instead of 2,254 km/s); J0927−6335 (2,522 km/s) is too fast for two ≤ 1.05 M☉ white dwarfs and needs a ≥ 0.95 M☉ donor with a 1.2 M☉ accretor (Bédard et al. radii; El-Badry et al.'s donor limits reproduced to 0.01 M☉ at 20,000 K). Orbits in McMillan (2017) with galpy: D6-2's path passes 0.95° from the centre of G70.0−21.5 92,000 years ago in all 600 draws (Shen et al.: 0.9°); by chance 3.2% for its track, 6.2% for any D⁶ star (20,000 longitude shuffles of Green's 310 remnants). Forward model of Gaia's selection with real parallax errors: bright, fast D⁶ survivors come from ≤ 0.8% (0.4–1.3%) of Type Ia supernovae (seen 2.2 Myr after launch); faint ones like D6-2 ≤ 7.8%, or 37% if they fade within D6-2's 92,000 years (El-Badry et al.: 3–5% and ~100%; Shen 2025: ~2%). Plus an offline 3D explorer flying all 17 stars out of the Milky Way. |
| M1 Visual binary orbits | **done** — [50-binary-orbits/](50-binary-orbits/). Bond et al.'s 2,367 measures of Sirius B (1862-2016) fitted from scratch (blind Thiele-Innes search), then with each of 117 visual observers' personal equation marginalised (MCMC): a = 7.4948 ± 0.0039″, e = 0.5914; Bond's Table 4 reproduced exactly with his weights. M_A = 2.057 ± 0.023, M_B = 1.015 ± 0.011 M☉ (Bond 2.063/1.018); the period is model-dependent (50.13-50.19 yr) but moves the masses only 0.3%. Sirius B is in Gaia DR3 (374.5 mas): its orbital curvature biases the 5-parameter parallax, and a scanning-law simulation conditioned on its proper motion corrects it to 379.0 ± 2.4 mas (Hipparcos 379.2; Hipparcos's own proper motion carried a −19.5 mas/yr leak from Gatewood's 1978 orbit). Photometric radius 0.0079 R☉, on the C/O mass-radius relation; GR redshift 81.8 ± 1.0 km/s vs measured 80.65. Mass-luminosity relation from 205 ORB6 grade 1-2 pairs (Sun-like star 1.03 ± 0.02 M☉; only 55 of 630 orbits have a clean Gaia parallax). Gallery animation + interactive page of six famous orbits. Spec corrected (ORB6 https and column overflow, Sirius B in DR3, download ~40 MB + 570 MB optional). |
| M2 Binary census | **done** — [51-binary-census/](51-binary-census/). 32,487 Sun-like stars within 100 pc, 37.7% with some Gaia DR3 sign of a companion. Gaia's detection simulated star by star (gaiamock for astrometry; an RV model calibrated on SB9, for which DR3 published only 437 of the 755 orbits its documented cuts allow): raw counts in four channels match Moe & Di Stefano within 16% and Raghavan within ~25%; corrected, 0.10–0.14 companions per star below ~9 yr and 0.21 beyond ~150 yr. Raghavan's log-normal rebuilt from his Table 13 (5.03, 2.28). Circularization period 9.7 d (Gaia), Bashi's temperature trend reproduced, no thin/thick-disc age difference; DR3's double-lined eccentricities are unphysical below 3 d. Twin excess fades with separation (half El-Badry+2019's at < 350 au). Wide-binary fraction 4.9% (M dwarfs) → 15% (Sun-like). |
| M3 Wide binaries and MOND | **done** — [53-wide-binaries/](53-wide-binaries/). El-Badry et al.'s Table 1 counts exact; Chae's 26,615 exact (two cuts found only in his code); Pittordis & Sutherland's search rebuilt in 8.86 M Gaia DR3 stars, Banik et al.'s sample 8,228 (8,611), Cookson et al.'s 1,220 (1,421). Own QUMOND solution (η = 1.4228) and orbit libraries for nine force laws; forward model with Gaia-fitted companion wobble: α = 0.003 (+0.014 −0.012) with Hwang's eccentricity prior (Banik et al. 0.024), Newton over MOND 19σ; with γ free α = 0.10 ± 0.03 because γ runs to ≈ 3 (MOND's fit to the prior edge, unlike theirs), but the pairs' own v–r angles, modelled with companions, prefer γ ≈ 1.5–2, where α = 0.03 and Newton leads by 12σ. Chae's method reproduced (1.55; his 1.49 ± 0.07) and injection-tested on mock skies of his own pairs: Newton 1.30–1.52 depending on the fitted companion population, MOND 1.78–2.11, so with companions fitted to his own sample the sky lies midway and the statistic alone cannot decide; ideal Newtonian binaries read 1.14 (per-pair Hwang eccentricities +0.25, his scale-free reference −0.1). The forward model on his own sample: α = 0.09–0.18 (2–4% speed-up), MOND out at 6–8σ. His sample's median-ṽ rise vanishes within velocity-precision classes. 3D, Chae et al. 2026's 36 pairs: γ = 1.84 (+0.21 −0.18), 1.14 without the six pairs too fast to be bound; Saad & Ting's 1.00 not reproduced (their error model gives 1.76 in our likelihood). Cookson et al.'s Δχ² = 12.6 not reproduced on our rebuild. Spec corrected (Chae et al. 2026; eccentricity watch-out). Offline explorer. |
| M4 Triple stars | **done** — [52-triple-stars/](52-triple-stars/). Multiple Star Catalog (2026-06-19): 819 of 854 inner/outer pairs with both orbits and masses lie above Mardling & Aarseth's limit (the 35 below are catalogue problems); TIC 290061484's period ratio is 1.7× its limit. Tokovinin 2017's spin-alignment statistics reproduced exactly from his tables (C = 0.61/0.23/0.18 by primary mass; 80% of 54 two-orbit triples within 70°) and repeated on the 2026 catalogue's 107 (85% within 80°); 33 of 42 binaries below 2.9 d have tertiaries (Tokovinin et al. 2006). TIC 290061484: 201 eclipses timed in TESS FFIs (median 0.65 min); the 24.52 (24.33–24.73)-d outer period found from the timing wave alone; a four-body REBOUND fit gives B's conjunction period 24.492 d (Kostov et al. 24.489 ± 0.004), and Kostov et al.'s own solution fits our times to 3.3 min after a 0.12% code-to-code recalibration of its inner period (timing alone moves e_out to 0.236 and m_B to 8.24 M☉, but would make B miss the pair). Eclipses that switch off and on, fitted rather than only simulated: one precessing-orbit model on 1,063–3,391 Harvard plates each (DASCH DR7 API), every TESS sector, Hipparcos and ASAS-3, with N-body checks of 300 draws per solution. V907 Sco: J = 27°, 60.7-yr wobble, eclipses end 2044.5 (Zasche et al. 2023: 26°, 63.5 yr, 2049; their forecast trails TESS 2023–25 by 0.5–1°). HS Hya: swept edge-on in 1964, 530-yr wobble, m3 = 0.54 M☉, eclipses return ~2132; Vokrouhlický & Zasche's retrograde solution (713 yr, ~2199) at odds 1:35. SS Lac: 625-yr wobble (Torres 2001: 600) in a narrower cone (tilt rate 0.095 vs 0.144°/yr), eclipses faded out ~2024 and return ~2092, not after 2200. Single nested-sampling runs get trapped (six HS Hya runs found four different geometries); SS Lac's free fit sat in the Kozai window and every N-body draw collided, while a Kozai-quiet fit found a better solution. Algol: point masses collide after 851 yr, GR no help, tidal + rotational bulges quench Kozai–Lidov (e < 0.004). Triples pushed past the limit break up as the m⁻³ law says (33/46/21%); MEGNO maps find Kozai chaos far beyond the limit near 90°. Proxima: P = 553 (+74/−55) kyr, e = 0.49 (Kervella et al. 547 kyr, 0.50); Gaia DR3 + Akeson et al.: 495 kyr. Bugs found and fixed: WHFast safe_mode=0 + exact_finish_time drift (46 min of eclipse timing by 2022), a 246-day finite-difference step on Julian dates, escaper mis-identification in a third of break-ups. Keeper animation and an offline explorer. |
| M5 Contact binaries and mergers | **done** — [54-contact-binaries/](54-contact-binaries/). Own Roche ray tracer (lobe volumes within 0.8% of Eggleton; with Molnar et al.'s parameters KIC 9832227's model amplitude is 0.194 mag vs Kepler's 0.195). ZTF's 369,578 contact binaries stop at P = 0.2272 d, the same in every apparent-magnitude slice (0.2248–0.2289 d); Gaia DR3's contact-like EBs give 0.2247 d, but Gaia's period search stops at 0.2 d; within 1 kpc the distribution peaks at 0.28 d (Rucinski: 0.27); the sub-edge systems are red-dwarf pairs (median M_G 8.2). Chen et al. 2018's W1 PL relation reproduced exactly (unweighted OLS) from their table; with Gaia DR3 parallaxes the same stars give −6.13 ± 0.17 log P − 0.19 (scatter 0.163 → 0.139 mag, 6.4% distances), predicting 41,386 independent ZTF parallaxes to 1.0003 ± 0.0001; ZTF's own slope depends on distance (saturation, parallax cuts); 261,725 contact binaries mapped with a Wesenheit relation. KIC 9832227: 49 own timings 1999–2025 (NSVS, SuperWASP, Kepler, ASAS-SN, AAVSO, TESS; 1.98 min rms vs Kovács et al. 2019); the NSVS archive labels MJD as "JD − 2400000", and reading it literally moves the 1999 time by 60.55 min (half a day minus one orbit), with which an exponential merger law fitted to Molnar et al.'s data gives t0 = 2021.3 (+0.9/−0.7) (theirs 2022.2 ± 0.7); with the right time no merger date fits; the period is 1.05 s shorter in 2025 than in 2011 (Ṗ = −2.0e−9 plus a 17-yr wave), and Kovács et al.'s 2019 model misses 2024 by 55 min. V1309 Sco (OGLE photometry is public; spec corrected): our 17 periods match Tylenda et al.'s to 0.13 s; −1.23% in six years; t0 = 2010.19 ± 0.59 (theirs 2010.10); the rise began 2008.22; the eclipse signal faded (0.12 → 0.035 mag) as a one-hump signal grew. T CrB: 950,487 AAVSO observations through the public VSX API (no account; spec corrected), 1866 to 2026-10-09: not erupted (V 10.29, B 11.72); 1946 t2/t3 = 3.3/5.0 d (Schaefer 3.0/5.0); secondary maxima at day 159/148; AAVSO's pre-1946 record is essentially Leslie Peltier's; 0.37-mag ellipsoidal wave whose irradiated-face minimum is 0.072 ± 0.009 mag shallower; recurrence alone gives 17% for the rest of 2026 (Pei et al. 2026's 30.2% reproduced for their date); false "eruption" reports filtered by a same-day observer consensus; one misdated 1946 record flagged. Plus an offline build-a-contact-binary explorer. |
| M6 White dwarf pairs for LISA | **done** — [55-white-dwarf-pairs/](55-white-dwarf-pairs/). ZTF J1539+5027 timed from scratch with our own forced photometry on 2,641 ZTF difference images (2018-2026) and 121 public PTF images (2009-2012): Ṗ = (−2.3670 ± 0.0024) × 10⁻¹¹ s/s (ZTF alone −2.3687 ± 0.0048), 0.4σ from Burdge's 2024 HiPERCAM ephemeris (Shaw et al. 2024) and 1.1σ from Burdge et al. 2019; eclipse times within 0.4-0.8 s of theirs. Chirp mass 0.3183 M☉ if pure GR; the decay is 1.09 ± 0.08 × GR's prediction from the published masses (tides expected to add 7%). LISA: all 18 of Kupfer et al.'s verification binaries reach SNR 7 within 3 months against instrument noise (7 with Robson et al.'s early-mission foreground); fates from Marsh et al.'s stability criteria; an observationally driven Galaxy holds 28 million pairs in the band, ~79,000 resolvable. Interactive LISA explorer and an orbit animation. |
| N1 Orbital precession | **done** — [57-orbital-precession/](57-orbital-precession/). Our own REBOUND/REBOUNDx integrations from DE440's J2000 state (1550-2650, 2-day sampling, IAS15; 37 runs) split Mercury's in-plane perihelion rate force by force: planets 532.289″/cy (Venus 277.49, Jupiter 153.82, Earth-Moon 90.86, Saturn 7.22, Mars 2.48), GR 42.9807 (EIH, single-Sun 1PN and Nobili-Roxburgh potential agree to 8e-6), solar J2 +0.0278 and Lense-Thirring −0.0020 (both equal their formulas) → 575.2956 vs DE440 575.2969 measured the same way; the 0.0013 gap is the asteroid belt (DE440's 343 masses: +0.0015). Park et al. 2017's own fit (AD 1000-3000, quadratic + 14 periodic terms) on our six-millennium run gives 575.3094 vs their 575.3100 ± 0.0015 (their Q −0.044781 reproduced as −0.044790). Frame trap: in the fixed J2000 ecliptic ϖ advances 571.93″/cy, on the ICRF equator 560.93 (the dogleg (1 − cos i)Ω̇). Newton-only Mercury falls behind DE440/DE441 by 43.01″/cy (279″ and 77,000 km by 2650, 21.5′ by AD 5000); with GR + Sun within 0.07″ over 1000 BC-AD 5000, the rate itself sliding from 578.0 to 572.7″/cy. Per-planet shares compared four ways with Park et al.'s Table 3 (pairs, leave-one-out, Gauss's 1818 averaging to 0.003-0.35%, Laplace-Lagrange 2-5% high); Jupiter's 0.54″/cy pair-vs-Gauss gap shown to be second order in its mass. History 1859-2017 (Le Verrier 38″ to MESSENGER 42.9799 ± 0.0009; Newcomb 1895 is 41.24 ± 2.09, not 43). Precession ladder of 14 measured orbits, x = GM/c²a(1−e²) from 3.6e-10 (LAGEOS II) to 1.9e-4 (S2): Baroch et al. 2021's TESS binary rates re-derived from their LaTeX timing tables (A = (2+2γ−β)/3 = 1.000 ± 0.012 vs 1.002), DI Her 0.96 ± 0.14 with tilted spins (−0.18 if upright), pulsars weighed without ω̇ (Hulse-Taylor 1.0008 ± 0.0007, double pulsar 1.000005 ± 0.000025), S2 GRAVITY 1.10 ± 0.19; our 1PN integration of GRAVITY's orbit moves the 2026 apocentre 0.782 mas (GRAVITY: 0.78). Labelled simulation: rosette movie with GR ×10⁶ (Nobili-Roxburgh potential; REBOUNDx's iterative 1PN force fails to converge there). Spec corrected (frame advice, download ~125 MB). Full rebuild ~25 min (six-millennium EIH run 18 min), from cache 1 min offline. |
| N2 Light bending | **done** — [56-light-bending/](56-light-bending/). Every table of Dyson, Eddington & Davidson 1920 and of the Lick Observatory's Wallal bulletins (1923, 1928) transcribed from the scans: ~13,500 numbers, each read twice blind (0 digit disagreements), with the check scripts in transcribed/checks/. The 1920 reduction replicated exactly (1.980″), which proves two Table II misprints and a misprinted Sobral longitude. Joint least squares: Sobral 4-inch 1.99 ± 0.14″ (1.6σ from Einstein, 7.7σ from Newton), Príncipe 1.74 ± 0.39″, the 1920 choice 1.96 ± 0.14″; the astrographic 0.96 ± 0.23″ is 3.8σ from the 4-inch (its scale and deflection trade as α + 0.24e). Table I matches Gaia DR3 carried back to 1919 to 1.7″ rms; the paper's Plate 1 is dated by JPL's Moon to 11:59:50 UT ± 7 s, which makes it 4-inch plate III. Wallal: five-foot cameras 1.83 ± 0.22″, 15-foot 1.69 ± 0.20″, combined 1.75 ± 0.14″ (Campbell & Trumpler: 1.75 ± 0.09″ p.e.). Cassini's raw 2002 tracking data from PDS (1.3 GB; own BCRS Doppler model, X/X + Ka/Ka plasma-free combination, the navigation team's trajectory): γ − 1 = (−0.5 ± 2.6) × 10⁻⁴, against Bertotti et al.'s (2.1 ± 2.3) × 10⁻⁵ and JPL's closed-loop ± 5.2 × 10⁻⁵. A century of γ from N7's table; interactive plate lab (artifact https://claude.ai/artifact/FksLym6sN7y5LupavdQDTi). |
| N3 Gravitational redshift | **done** — [58-gravitational-redshift/](58-gravitational-redshift/). Ten years (2016-2026) of ESA's public 5-min orbits and clocks for the two Galileo satellites stranded in eccentric orbits: with the IGS term −2r·v/c² added back each hydrogen maser swings ±0.38 µs once per 12.9-h orbit exactly as GR predicts (convention verified on the data; own proper-time model with J2 and tides, checked to 0.05 ps on synthetic orbits). Global fit of Delva et al.'s α over 2,605 days of E14's maser and 3,382 of E18's, with Sun-locked nuisance signals (E18 carries a 160 ps once-per-orbit signal locked to the Sun, not the perigee) and errors calibrated by fake-perigee null tests: α = (−2.6 ± 4.1) × 10⁻⁵ (clock noise 1.2, model choice 1.3, orbit 3.2 from ESA's SLR summaries, magnetic 1.9 adopted), 3.4× tighter than Gravity Probe A and consistent with Delva (0.19 ± 2.48) and Herrmann (4.5 ± 3.1); more data did not beat their systematics. E14's raw 2016-17 value (+11 ± 6) matches their raw values and becomes −0.7 ± 4.1 with Sun-locked terms. The beyond-Kepler relativistic term the IGS convention leaves in (J2 orbit bending + J2 potential + tides, 110-150 ps) found in all 8 clocks at 1.02 ± 0.04. GPS: 3,852 satellite-days, GR's gain 38.574 µs/day vs drift-minus-factory-offset 38.556. Sun (IAG atlas, own core fits, NIST lab scale): 622 ± 26 m/s vs 633.1 (643 ± 13 on González Hernández et al.'s 15 lines, theirs 639 ± 14). Sirius B: M1's 81.8 ± 1.0 km/s vs Joyce's 80.65 ± 0.77. Redshift ladder of 23 rungs (literature values checked against the primary papers). E14's maser degrades from 2024 (clock history from the data). Animation of E14's clock. Spec corrected (data route, size ~3.4 GB streamed / 0.45 GB kept, new watch-outs). |
| N4 Geodetic precession and frame dragging | **done** — [61-gyroscope-precession/](61-gyroscope-precession/). **Frame dragging measured from public data**: 4.1 years (Aug 2022 – Oct 2026) of the ILRS analysis centres' weekly LAGEOS and LARES 2 orbits (EDC, anonymous FTP; 8 centres, 3,435 files) refitted with an own numba IERS-2010 force model (GOCO06s, solid/FES2004 ocean/pole tides, DE440, SRP with shadow, Schwarzschild + de Sitter; frame dragging off); weekly node jumps summed over the supplementary pair, J2 from their difference, periodic tide/SRP terms fitted: μ = 1.002 ± 0.018 (stat) ± 0.011 (syst) of GR's 61.35 mas/yr (47σ), seven centres 0.999–1.006 (NSGF's week-to-week frame rotations give 0.55 ± 0.54), injection-test gain 0.997; Ciufolini et al. 2026: 1.001 ± 0.002. Trap found: fitting C_R weekly absorbs 0–19% of frame dragging (synthetic) and fakes −21 mas/yr inclination drifts; C_R is held fixed. De Sitter (17.2/17.0 mas/yr) degenerate with LT: switched off, μ = 1.574 (expected 1.560). By-product: J2 from the pair tracks CSR's monthly C20 (r = 0.88) with a 14σ offset of 1.14 × 10⁻¹⁰ traced to the IERS FES2004 file's 18.6-yr equilibrium ocean tide being anti-phased to the solid tide; re-run in phase: −0.05 ± 0.07. GP-B recomputed along its orbit: frame dragging −39.1 (−39.2), Sun +7.4/−16.4 (+7.3/−16.2), geodetic −6,597 to −6,620 depending on mean radius (2.35 mas/yr per km; −6,606.1 ↔ 7,023.6 km), lunisolar cancellation of the J2 node drift (0.30° → −0.025°). Pulsar spin-precession rates recomputed (5.074, 2.234, 0.514 °/yr); B1913+16 and J1141−6545 have no measured rate. TLEs 265–340 m off cross-track. Spec corrected (data route, ~0.16/0.9 GB, pulsar list, watch-outs). Full campaign ~8,700 arc fits (~4 h on 8 cores), run.py 4 min, offline-verified. |
| N5 Falling into a black hole | **done** — [62-falling-into-a-black-hole/](62-falling-into-a-black-hole/). Labelled simulation: a 76-s movie of a radial fall into Sgr A* (H1's 3.96 million M☉) from 1,000 r_s to 0.01 r_s, the last 26 s (horizon to singularity) in real time, against the real sky: 1.25 million Gaia DR3 stars (G < 11, Hipparcos below V = 5) lensed one by one plus the summed light of all 1.8 billion fainter ones, aggregated on the archive's server into 3.4′ HEALPix pixels. Exact ray map from Carlson elliptic integrals in the faller's own (Gullstrand-Painlevé) frame, inside the horizon too, checked against 40-digit quadrature and brute-force geodesics built on SymPy-derived Christoffel symbols (≤ 10⁻⁸ rad), the weak- and strong-deflection series (Bozza 2002) and Synge's shadow formula. Colours: blackbodies shifted T → gT. Falling, the dark region is 42.1° in radius at the horizon (hovering observers see it fill the sky); the sky ahead is blueshifted ×3.875 at the shadow's edge, ×½ straight behind; near r = 0 the outside universe becomes a thin blue band round the waist (as Hamilton describes). Clocks: 13.3 min from 10 r_s to the horizon, 26.0 s inside (π GM/c³ = 61 s at most); a far-away friend sees the redshift grow ×e every 4.001 GM/c³ (78 s) and the watch freeze. Tides: whatever the mass, a 2 m body feels 10 g 0.095 s before the singularity ((2/3)√(L/a)); you cross the horizon intact only above ~14,500 M☉ (Cyg X-1 tears you apart 78 r_s out). Kerr: frame dragging (GW250114's remnant a = 0.68: 0.4 turns, horizon 101 Hz; Cyg X-1 a = 0.9985: 5.8 turns), ISCO plunges, the ergosphere. SymPy: GP = Schwarzschild, flat slices, Kretschmann 48M²/r⁶, tidal tensor. Also Penrose and river diagrams, clocks and tides for M87*, Sgr A*, GW250114, Cyg X-1 and Gaia BH1. 17 tests pass; runs offline once cached. Spec corrected: data route (Gaia server-side flux map, not a hips2fits image; I7 not built, so GW250114 values from LVK 2025), ~0.36 GB downloaded (0.22 GB kept), watch-outs. |
| N6 Killing vectors | **done** — [63-killing-vectors/](63-killing-vectors/). Killing vectors counted exactly by Taylor-expanding Killing's equation about one point in rational arithmetic (a rigorous upper bound), matched by explicit polynomial solutions verified by substitution: flat 10, Schwarzschild 4, Kerr 2, FLRW 6 for generic/matter/radiation/Planck ΛCDM a(t) (none timelike; a(t)∂t conformal), de Sitter 10, Gödel 5, Einstein static 7, plane wave 5. Kerr's hidden symmetry found by the same count (5 invariant Killing tensors, one irreducible; Zipoy–Voorhees 4) and verified as the Penrose–Floyd Killing–Yano square. ISCO, photon sphere, √27, 5.719-42.265% efficiencies from the Killing-reduced radial equation (BPT to 1e-13). S2 (H1 data and chain): areal velocity in 13 independent seasons χ² = 10.7/13 (2,039 ± 32 vs 2,025 mas²/yr) while speed changes 13-fold; energy in 12 seasons χ² = 9.9/12; precession from E, L 11.70' ± 0.07' (GRAVITY 2020-2024 all consistent). Kerr shadow 9.66-10.39 GM/c²; Sgr A* 51.2 ± 0.3 µas vs I5's 53.9 ± 1.3 ring. Hawking T, entropy, lifetimes for 32 atlas black holes + 746k quasars; GW250114 area 236,500 → 374,200 km² (LVK samples). D7's β reproduced exactly; FLRW photon a·E constant to 6e-8. Zipoy–Voorhees orbits looked regular (KAM; documented negative result). Kerr orbit explorer page. 34 self-tests. Spec corrected (data 27 MB, I7/I4 substitutions, S2 test design). |
| N7 Tests of GR scorecard | **done** — [64-gravity-scorecard/](64-gravity-scorecard/). 46 tests on one footing (1σ fractional deviation from GR, or a violation parameter GR sets to zero, against GM/rc² from 3 × 10⁻¹² to 1/3): 35 published landmarks plus 11 atlas measurements (N3 ×4, I5 ×2, H1, M3, D10, N7 ×2). 36 lie within 1σ and 45 within 2σ (31 and 44 expected); the only outlier is the growth index, +3.4σ (Nguyen+23, built on Planck PR3 lensing). 625 published values transcribed with each paper's exact sentence and confidence convention into literature/checks/ (182 confirmed the first draft, 19 corrected it, 388 added, 36 secondary-only). Errors found in Will (2014): Lebach+95 quoted at twice its error; J1738+0333's α₁ 95% interval turned into a 7 × 10⁻⁵ "limit" (the interval's full width); early eclipses plotted with probable errors. Our own lunar laser ranging: JPL's published range model (Williams & Boggs 2020: IAU 2006/2000A, EOP 20 C04 + ORTHO_EOP, solid, pole and lunar tides, Mendes–Pavlis + FCULa, Sun and Earth Shapiro delays, TDB scale terms; exact-time Chebyshev evaluation) fits all 3,846 APOLLO normal points (2006-2020, 3.1 mm median precision) to 2.1 cm rms with EPM2021. Residual consistency tests: γ − 1 = (−2.7 ± 5.6) × 10⁻³, Nordtvedt η = (1.7 ± 1.4) × 10⁻⁴, and APOLLO's plate motion matches GPS station P027 to within 1.2σ. Unresolved and documented: with DE440 a 7-cm anomalistic-month residual remains, the Moon 1-2 m (~1.3 ms) ahead along its orbit compared with EPM2021, DE421, DE430 and what APOLLO prefers year by year (DE441 is ahead too). Simulations: Nordtvedt effect 12.98 m × η (REBOUND; literature 12.9-13.3); triple pulsar 34 ns per Δ = 10⁻⁶ (Archibald+18: 30), white-noise limit 1.35 × 10⁻⁶ (95%) from Voisin+25's 12,474 public TOAs vs their 1.46-2.29 × 10⁻⁶. Also: WEP timeline Newton → MICROSCOPE, γ timeline 1919 → Cassini, γ–β plane, mass–mass diagrams (2PN ω̇), GR prediction table, offline interactive scorecard page. 18 self-tests; runs offline once cached; a forced rerun reproduces every number. Spec corrected (Voisin 2025, the LLR inputs, Will 2014 errors, residual-test caveat). |
| O1 Friedmann universes | **done** — [65-friedmann-universes/](65-friedmann-universes/). Exact solutions of the Friedmann equation for every kind of universe (turning points from the roots of P(a) = a⁴E², time integrals with the square-root singularities removed by substitution; checked against every closed form, a direct integration of the acceleration equation to 10⁻¹⁰, and Carroll, Press & Turner's boundaries to 10⁻⁹): a gallery of 12 universes, a map of fates with a(t) glyphs, and an interactive offline build-a-universe page. Curvature free: Pantheon+ alone reproduces D1 (Ωm = 0.294 ± 0.054, ΩΛ = 0.609 ± 0.081); DESI DR2 BAO alone reproduces DESI's own chain to 0.01σ (Ωk = +0.025 ± 0.040; χ² at DESI's best fit 9.959343 vs 9.9593429); the compressed Planck PR4 CMB alone is a thin geometric-degeneracy band from (Ωm 0.145, H0 99) through flat at Ωm 0.3125, H0 67.4 to (1.2, 34). BAO + CMB: Ωk = +0.0021 ± 0.0012, H0 = 68.61 ± 0.33 (DESI DR2 with the full CMB: +0.0023 ± 0.0011, 68.50 ± 0.33; an independent MCMC agrees with our grid + Laplace marginalisation to 0.02σ); all three: Ωk = +0.0023 ± 0.0012, so a curved universe has a radius above 209 billion light-years. Age 13.71 ± 0.05 Gyr curvature free, 13.787 ± 0.016 flat (Planck + BAO 13.787 ± 0.020); with the Cepheid H0 the same universe is 12.77 Gyr old, 1.0σ under B6's oldest stars, while Einstein–de Sitter is 8.87 Gyr, 5.7σ under (the 1990s age crisis). Epochs: matter–radiation equality z = 3383 (51.6 kyr), CMB z = 1090 (372 kyr), acceleration from z = 0.66 (6.2 Gyr ago), dark energy overtakes matter at z = 0.32. New result: the CMB's own photons cap any bouncing (no-Big-Bang) universe at z ≤ 9.34, so quasars at z = 7.6 need Ωm < 0.0017 and the z = 14.4 galaxy and the CMB are impossible: our universe had a Big Bang. Found that CAMB's physical_time() is accurate only to ~2×10⁻⁵ (0.28 Myr). Spec corrected (download ~0.14 GB, Planck age label: 13.787 is +BAO; B10 not yet built, so only B6's ages are used). |
| O2 Light cone and horizons | **done** — [66-cosmic-light-cone/](66-cosmic-light-cone/). Own Friedmann integrator (photons, neutrinos incl. one 0.06 eV, matter, Λ or w0–wa), within 2 × 10⁻⁶ of astropy and 1.5 × 10⁻⁵ of CAMB. Planck 2018: particle horizon 46.20 ± 0.07 billion ly, event horizon 16.58 ± 0.15, recession = c today at z = 1.476, smallest angular size at z = 1.592, event horizon at z = 1.852, 4.62% of the observable universe still reachable (97% of it recedes faster than light). Davis & Lineweaver 2004's Fig. 1 read back from the vector EPS in their arXiv source: worldlines match to 0.001 billion ly, light cone/Hubble sphere/event horizon to ≤ 0.04; drawn without radiation, conformal times 0.05–0.16 short; their text numbers (z 1.46 and 1.8, CMB 3.2c and 58.1c) reproduced. 1.4 million real objects placed on the light cone (SDSS galaxies and DR16Q quasars, Pantheon+, 282 GWOSC mergers, 3,223 JADES DR4 galaxies, MoM-z14, the CMB): 61% of DR16Q quasars recede faster than light today and 42% are beyond our event horizon; the CMB's gas emitted from 42 million ly and will only ever be seen up to age 0.77 Gyr. Measured: DESI's Lyα forest receded at 1.356 ± 0.029c at emission (12σ; ruler- and H0-free), DESI alone puts v = c at emission at z = 1.60 (+0.08/−0.07) and sees D_A falling at z = 2.33 (12σ); Pantheon+'s 7 supernovae beyond z = 1.48 recede at 1.114 ± 0.043c today. Finding: D&L's 23σ against the special-relativistic reading of redshift comes from curves drawn with D_L = D, not their text's D(1+z) (reproduced: 22.7σ); Pantheon+ rejects the text recipe at 22σ and the self-consistent SR (empty Milne) universe at 7σ. Offline interactive explorer (follow a photon; observe from the future). Spec corrected (data route, ~5 MB, D&L watch-outs). |
| O3 Olbers' paradox | not started |
| O4 Recombination | not started |
| O5 21-cm cosmic dawn | not started |
| O6 Cosmic dipole | not started |
| P1 Inflation | **done** — [67-inflation/](67-inflation/). Eight inflaton potentials solved exactly (background + Mukhanov–Sasaki modes; power-law inflation's exact n_s, r reproduced to 10⁻⁸; first-order slow roll overestimates n_s by up to 1.9 × 10⁻³, second order agrees to 2 × 10⁻⁵), with N* from reheating solved self-consistently (Starobinsky instant: 55.7, matching Planck 2018 X eq. 47). Our own re-implementation of the BICEP/Keck BK18 likelihood (Hamimeche–Lewis, 11 maps, 594 bandpowers) matches Cobaya's χ² to 0.008 rms over 120 samples and gives r < 0.0361 (BK18: 0.036), A_d = 4.48 +0.80/−0.71, β_d = 1.50 ± 0.13; with Planck 2018 r < 0.0337, n_s = 0.9654 ± 0.0040 (ACT's joint run: 0.0337). All 16 chains (Planck PLA grid and SPT-3G read by byte range from 11.8 and 6.2 GB zips; ACT DR6 from LAMBDA) reproduce their papers within 0.15σ. Verdicts (best point in N* = 50–60): φ, φ², φ⁴ beyond 3.3σ everywhere; natural ≥ 3σ; φ^(2/3) 3.2σ (Planck + BK18) but 2.0σ with ACT; Starobinsky 0.3σ / 1.7σ / 1.5σ (Planck + BK18 / P-ACT-LB + BK18 / CMB-SPA + DESI), 2.3σ with instant reheating against ACT and SPT; reaching ACT's n_s needs N* ≈ 75 (1σ: 64–67, only a stiff w → 1 phase). Horizon problem on the Planck sky: θ_H = 1.16°, ~39,000 patches; flatness tuned to 5 × 10⁻⁶³ at the Planck temperature without inflation; 61.1 e-folds needed. Offline interactive explorer. Spec corrected (SPT wording, λφ⁴'s r, data routes, BK-XX maps not public, watch-outs). |
| P2 Thermal history | **done** — [68-thermal-history/](68-thermal-history/). g*(T) from 10 TeV to 10 keV built from the 2026 PDG list. Ideal gases give 106.75 → 10.75 → 3.363/3.909 exactly. Our interacting build (a gas of 454 PDG hadrons → Borsányi et al. 2016's lattice QCD → perturbative QCD to g⁶, whose unknown constant matches the lattice at 1 GeV for q_c = −3170 against their −3000) gives 104.8 at 10 TeV and 3.3828/3.9302 today. It lies inside Saikawa & Shirai's band at 82% of temperatures and within 0.8σ above 10 MeV; the ideal gas is 2.1× too high just above the 156.5 MeV crossover, and pions alone hold a fifth of lattice QCD's strongly interacting matter at 150 MeV. Neutrino decoupling solved with Escudero's equations, QED corrections and the 2025 electron-mass rates: N_eff = 3.0443 (Bennett et al. 2021: 3.0440 ± 0.0002), reproducing Escudero 2020 and Escudero et al. 2025 to every printed digit. Calendar: electroweak crossover 9.4 ps, quarks bind 19.5 ± 0.6 µs, neutrinos decouple 0.12–0.34 s, n/p freezes 1.4 s, deuterium survives 106 s, equality 51,300 yr (z = 3392), CMB 372,300 yr, age 13.787 Gyr (clock matches Escudero's to 0.01%). QCD softening w = 0.229 at 178 MeV (Saikawa & Shirai 0.229 at 176 MeV), with 2.6 M☉ inside the horizon; horizon entry then maps to 3.4 nHz today (pulsar-timing band). Spec corrections above: ~3 MB, interacting check values, horizon mass, watch-outs. Plus the Cosmic Thermometer page (offline, and a private artifact). |
| P3 Cosmic neutrino background | **done** — [69-neutrino-background/](69-neutrino-background/). Relic neutrinos today: 336.0 per cm³ at 1.945 K (339.8 with the heat caught at decoupling); a 50 meV state drifts at a median 2,850 km/s, too fast for the Galaxy to gather. Decoupling re-solved in independent code (Escudero's equations, QED to O(e³), electron-mass rates): N_eff = 3.0443, every row of Escudero 2020's Table 1 reproduced, identical to P2's to five decimals (Bennett et al. 2021: 3.0440). CAMB 2.0.4 fits with our own numpy likelihoods (Planck plik_lite, low-ℓ TT and EE, Planck 2018 and ACT DR6 lensing, BAO 2018, DESI DR2 BAO), each reproducing the published χ² at the published best fits (low-ℓ TT 23.2515 and EE 396.0499 exactly). N_eff by profile likelihood: D3's TT 2.62 ± 0.24 (too low without ℓ < 30, where n_s slides to 0.94), Planck TT + low-ℓ 2.94 −0.27/+0.28 (Planck 3.00 ± 0.28), TT,TE,EE 2.91 ± 0.19 (2.92 ± 0.19), + lensing + BAO 2.97 ± 0.17 (2.99 ± 0.17), Planck + ACT lensing + DESI DR2 3.14 ± 0.16 (Elbers et al. 3.23 ± 0.18 with CamSpec); N_eff ≈ 4 excluded at > 5σ, and SH0ES's H₀ would need N_eff = 3.9 (4.5σ). Σm_ν (95%): D3's TT < 0.36 eV, + DESI DR2 < 0.109; Planck rungs < 0.547/0.268/0.249/0.119 (Planck 0.54/0.26/0.24/0.12); Planck + ACT lensing + DESI DR2 < 0.0673 (DESI with the same Planck likelihood: 0.0691); w0wa + DES-SN5YR < 0.127 (0.129). The inverted ordering's minimum is disfavoured by Δχ² = 4.5 relative to the normal one's; each ordering on its own gives < 0.106/0.139 eV, so the lightest mass < 25/27 meV (Elbers 23/24); a parabola continued to negative mass puts the best fit at −0.038 ± 0.045 eV (Elbers −0.036 ± 0.043), 2.1σ from the normal minimum (1.1σ with w0wa). Answer keys read directly from DESI DR2's and Planck's chains (0.0642/0.0691/0.1291 and 0.1205 eV). KATRIN's m_β < 0.45 eV allows Σm_ν < 1.35 eV. Matter power at k = 1 h/Mpc: −7.5 to −8.5 f_ν linear, ~−10 f_ν non-linear. 50 self-checks; runs offline from cache. Spec corrected (data route, ~0.95 GB, the TT-only watch-out, DESI DR1's 0.082, the lightest-mass trap). Plus the Weigh a Neutrino page (offline, and a private artifact). |
| P4 WIMP freeze-out | **done** — [70-wimp-freeze-out/](70-wimp-freeze-out/). Own Boltzmann solver (Radau + an L-stable implicit stepper) with Saikawa & Shirai's g*(T): thermal ⟨σv⟩ = 2.11 × 10⁻²⁶ cm³/s at 100 GeV (2.11–2.39 for 10 GeV–10 TeV, peak 4.6 at 0.35 GeV), reproducing all eight of Saikawa & Shirai's 2020 tables to 0.09–0.13% and Steigman et al. 2012 in their own set-up to <1% outside 3–5 GeV; SDB12's 2.2 → 2.1 is mostly Ω h² 0.11 → 0.120. Freeze-out at x_f = 23.3 (4.3 GeV, 15 ns); Kolb & Turner's formula within 3% except in the QCD era (×0.82). Neutrons with the same physics: weak freeze-out at 0.67 MeV (1.7 s), helium forms at 72 keV (246 s); Born rates + 8 reactions give Y_p = 0.2426, D/H 2.43 × 10⁻⁵, matching PRyMordial's equivalent mode to 10⁻⁵ (its full corrections: +0.0044 → 0.2469; observed 0.2453 ± 0.0034). Electrons freeze out at z ≈ 407 (x_e = 3.40 × 10⁻⁴ at z = 200; CAMB 3.37). Own Fermi-LAT limits from McDaniel et al.'s 14-yr per-dwarf likelihoods (their profiles reproduced to ΔTS < 0.05): thermal WIMPs excluded below 78 GeV (b b̄) and 66 GeV (τ⁺τ⁻) for 42 dwarfs; Ackermann 2015's 6-yr curves (read from the vector PDF) cross at 94/74 GeV. Higgs-portal singlet relic line: LZ 2024 leaves only 59.5–62.5 GeV between 20 GeV and 30 TeV. Spec corrections above (2.2 is for Ω h² = 0.11, ~30 MB, data routes, watch-outs). Plus the WIMP Freeze-out Lab page (offline; solver agrees with Python to 2 × 10⁻⁴). |
| P5 Antimatter | **done** — [71-antimatter/](71-antimatter/). **η two ways.** We ran PRyMordial ourselves (both rate compilations, 60-draw Monte Carlos of the key rates) on the PDG's 12 quasar D/H systems (2.508 ± 0.029 × 10⁻⁵; the PDG and Cooke et al. means reproduced exactly). η = (6.03 ± 0.06) × 10⁻¹⁰ with PRIMAT rates (Pitrou et al.'s Ω_b h² = 0.02195 reproduced at 0.02194) and (6.11 ± 0.16) × 10⁻¹⁰ with NACRE-II, against the CMB's 6.08 (D3) to 6.16 (ACT): 0.1–1.9σ, the spread set by the d + d rates. **In numbers:** 410.7 photons per cm³ against 0.2515 baryons per m³ (1.63 billion photons per baryon); above 1 GeV, one extra quark per 273 million quark–antiquark pairs (the spec's "billion" corrected). A symmetric universe would freeze out at 21 MeV with η = 3 × 10⁻¹⁹, two billion times too little matter. **Domains.** Cohen, De Rújula & Glashow 1998 reproduced with our own 1D fluid solution of their equations, their curves read exactly from the arXiv PostScript: J within 2–31% and T within 0–17% for 1 + z ≥ 30, spectrum within 3–9%. Redone for Planck 2018 with 3.93 photons per annihilation (Crystal Barrel channels; the Backenstoss spectrum matched to 1–4%) against COMPTEL's MeV background: domains must exceed d₀ = 178 Gpc. The other cases give 346 Gpc with full heating, 11 with annihilation ending at z = 49 and 3.5 with no heating; CDG's own curve against the final COMPTEL data gives 9.5. A single boundary within 11.8 Gpc of us would outshine the MeV sky. **Cosmic rays.** One leaky-box grammage fitted to AMS-02 B/C (7.60 g/cm² × (R/10 GV)^−0.503) predicts AMS-02's p̄/p at 0.96–0.99 above 20 GV with no free parameter; antihelium < 6.9 × 10⁻⁸ (BESS). **Bullet Cluster** antimatter fraction < 5.6 × 10⁻⁷ with Fermi-LAT (Steigman 2008 with EGRET: 3 × 10⁻⁶). **Sakharov.** J = 3.158 × 10⁻⁵ from the PDG 2026 fit; Jarlskog's measure is 1.2 × 10⁻¹⁹ at 100 GeV down to 7 × 10⁻²⁴, 9–13 orders short of Y_B = 8.7 × 10⁻¹¹. Sphalerons freeze out at 131.5 GeV (lattice paper: 131.7); v/T = 0.17 at 125 GeV, and the lattice finds a crossover. **Found:** Korsmeier et al. 2018's eq. 14 prints (10 − √s)², but only the first power reproduces their own Fig. 4; the square fakes a 2.5× p̄/p bump at 9 GV. CDG's Fig. 4 lower curve is the upper one shifted 1.58 decades, not 1.70. Spec corrected: quark ratio, COMPTEL needed for the domain test, PRyMordial standing in for K1, Jarlskog range, ~3 MB. Plus an offline domain explorer (also a private artifact). |
| P6 Quark–gluon plasma | **done** — [72-quark-gluon-plasma/](72-quark-gluon-plasma/). Our own hadron resonance gas (459 PDG 2025 hadrons and nuclei from Thermal-FIST's tables, quantum statistics by Gauss–Laguerre quadrature equal to the Bessel series to 10⁻⁹, Breit–Wigner widths, feed-down by ALICE's and STAR's conventions, charm-decay pions and kaons) fitted to 21 ALICE yields from π to anti-⁴He (Pb–Pb 2.76 TeV, 0–10%): **T = 155.7 ± 1.4 (fit) ± 0.7 (model) MeV** (1.81 × 10¹² K), μ_B = 1.9 ± 4.2 MeV, V = 5,440 ± 470 fm³ with Andronic's excluded volume; χ² 40/18, protons −2.8σ (the proton anomaly), K* −3.1σ (rescattering, not fitted). Andronic et al. 2018 (156.5 ± 1.5, 0.7 ± 3.8, 5,280 ± 410) reproduced on their own data (155.9, 0.8, 5,470; their Fig. 1 read back from the vector PDF): at their best fit our Ω and nuclei equal their model to < 1%, Λ and φ differ by 15%/9% through feed-down. Nuclei alone 161.4 ± 4.3 (Andronic 159 ± 5), p+d+³He 157.6 ± 2.9 (ALICE 156 ± 4). STAR's GCE fits reproduced (7.7 GeV 142.9/398 vs 143.8/400; 39 GeV 157.9/104 vs 159.9/105, THERMUS yields within 2%); own fits at 39–7.7 GeV trace freeze-out from (106, 159) to (410, 144) MeV on the lattice crossover line. T independent of centrality (155.4–159.3). Hadron gas = lattice EoS (HotQCD 2014, WB 2014) below T_pc, then overshoots; χ₂^B parts from the lattice just above T_pc. The universe passed our T at 19.8 ± 0.5 µs (P2's calendar; own Friedmann integration 19.9) with μ_B/T = 1.1 × 10⁻⁸ (≈ 2 eV); antibaryons gone at 36 MeV (0.52 ms). Direct photons: T_eff = 304 ± 11 ± 45 MeV (ALICE 304 ± 11 unsubtracted), the universe at 3.7 µs. 18 self-checks; offline. Spec corrected (photon slope wording, data routes, excluded-volume/charm/feed-down watch-outs). Plus the Freeze-out Thermometer page (the gas runs in the browser; offline, and a private artifact). |
| P7 Higgs in CMS open data | **done** — [73-higgs-boson/](73-higgs-boson/). 15 GB of 2012 CMS open data, verified against CERN's checksums. Dimuon spectrum: 24.1 M pairs (cut flow identical to ROOT's df102), nine particles from η to Z; φ, J/ψ, ψ(2S), ϒ(1S–3S) come out 0.9–1.3‰ light and Z → μμ 0.16% light against the simulation, electrons 0.4–1.0% heavy (the files are uncalibrated); A_FB changes sign at the Z (photon–Z interference). The CMS outreach Higgs analysis (record 12360) reproduced in all 9 plots × 36 bins, read back from its vector PDFs, including its `int` pairing quirk; found a double count in it: 54 of the 94 2e2μ collisions are in both data streams (3 of them in the two peak bins), so the published histogram's 4.3σ is 3.1σ with each collision counted once. Our selection (any number of leptons, lepton momenta calibrated on the Z peaks, misidentified-lepton background from same-sign data): 12 events at 121.5–130.5 GeV (3.7 ZZ + 6.7 Higgs expected); m_H = 125.08 +1.10 −2.77 (stat) ± 0.12 (syst) GeV (PDG 2026: 125.13), μ = 0.81 +0.48 −0.39, local 2.5σ at the known mass (20,000 toys: 2.51σ; a Standard Model Higgs gives 3.0σ on average), global 1.1σ over 110–160 GeV; on-shell ZZ 170 events vs 186 simulated. Theory: our matching and three-loop running reproduce Buttazzo et al.'s tables (λ(M_Pl) −0.01435 vs −0.0143), and SMDR (built from source) agrees on the running to 10⁻⁵. PDG 2026: λ < 0 above 4.9 × 10¹⁰ GeV (one-loop running: 2.4 × 10⁸); absolute stability needs M_t < 171.05 GeV, 3.3σ below the measured 172.60; leading-order probability of decay in our past light cone 10⁻⁸⁵², (Γ/V)^−¼ = 10²²³ yr (Andreassen et al. 10¹⁶¹ at M_t = 173.1). The Higgs field switched on at 159.6 GeV, 9.4 ps after the Big Bang (lattice data + P2). Also found: a g_Y typo in Buttazzo et al. and their g₃(M_t) 0.0018 high, their quoted 171.53 GeV bound inconsistent with their own formula (171.12), and a past-light-cone volume of 0.123/H₀⁴ (not 0.15). Offline vacuum-stability explorer (three-loop RGEs in the browser). Spec corrected (PDG 2026, ~15 GB, efficiencies, duplicates, calibration). |
