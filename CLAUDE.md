# Exohuntr - Claude Code Instructions

## Project Overview

Exohuntr is a Rust + Python exoplanet transit detection pipeline that downloads NASA TESS light curves and runs BLS (Box-fitting Least Squares) transit detection to find planet candidates. The Rust engine handles compute-heavy BLS, validation, and cross-matching in parallel via Rayon; Python handles data access, deep analysis, and visualization.

**Repo:** humancto/exohuntr
**GitHub Pages:** humancto.github.io/exohuntr

## Architecture

```
exoplanet-hunter/
├── Cargo.toml                          # Rust project config
├── src/
│   ├── lib.rs                          # Rust library root
│   ├── bls.rs                          # BLS algorithm, SNR, phase math
│   ├── validate.rs                     # 5 false-positive tests + scoring
│   ├── crossmatch.rs                   # Hash-based catalog cross-matching
│   ├── io.rs                           # CSV parsing, file discovery
│   └── main.rs                         # CLI: search, validate, crossmatch subcommands
├── python/
│   ├── download_lightcurves.py         # Fetch data from MAST/NASA
│   ├── analyze_candidates.py           # Phase-fold plots, cross-matching, reports
│   └── validate_candidates.py          # 5-test false positive validation pipeline
├── data/
│   └── lightcurves/                    # Downloaded CSV light curves
├── results/
│   ├── plots/                          # Phase-folded light curve PNGs
│   ├── REPORT.md                       # Analysis report
│   ├── VALIDATION_REPORT.md            # Scored validation results
│   ├── validation_results.json         # Detailed per-candidate validation
│   └── crossmatch_results.csv          # Known exoplanet cross-match
├── docs/                               # GitHub Pages site
│   ├── index.html                      # Interactive results page
│   ├── candidates.json                 # Detection data for web UI
│   └── *.png                           # Phase-fold plot images
├── tests/
│   ├── conftest.py                     # Shared Python test fixtures
│   ├── test_validate_candidates.py     # Python validation tests
│   └── test_analyze_candidates.py      # Python analysis tests
├── candidates.json                     # BLS output (intermediate)
├── Makefile                            # Build & run automation + `make test`
├── pytest.ini                          # Python test configuration
├── scripts/setup.sh                    # One-command setup
└── CLAUDE.md                           # You are here
```

## Testing

Run all tests (98 total: 66 Rust + 32 Python):

```bash
make test          # Run everything
cargo test         # Rust only (66 tests: bls, validate, crossmatch, io)
python3.11 -m pytest tests/ -v  # Python only (32 tests)
```

### Rust test modules:
- `bls::tests` — BLS period recovery, SNR estimation, phase math, median (18 tests)
- `validate::tests` — All 5 false-positive tests, scoring, integration (23 tests)
- `crossmatch::tests` — Catalog indexing, lookup, CSV loading (13 tests)
- `io::tests` — CSV parsing, NaN handling, file discovery (8 tests)

### Python test files:
- `tests/test_validate_candidates.py` — Validation functions, scoring (24 tests)
- `tests/test_analyze_candidates.py` — Plotting, cross-matching, binning (8 tests)

## Pipeline Steps

### Step 1: Download light curves

```bash
# IMPORTANT: Use python3.11, not python3 (which may be 3.9 without deps)
python3.11 python/download_lightcurves.py --candidates-only --limit 500 --catalog
```

- `--candidates-only` downloads unconfirmed TOIs from ExoFOP (best for discovery)
- `--mission tess --sector N` downloads by sector (note: sector search by number can fail in lightkurve; prefer --candidates-only)
- `--catalog` saves stellar parameters alongside light curves

### Step 2: Build and run Rust BLS engine

```bash
cargo build --release
# Using subcommand (recommended):
./target/release/hunt search -i data/lightcurves -o candidates.json --snr-threshold 6.0
# Backward-compatible top-level flags also work:
./target/release/hunt -i data/lightcurves -o candidates.json --snr-threshold 6.0
```

Aggressive search:

```bash
./target/release/hunt -i data/lightcurves -o candidates.json \
  --snr-threshold 5.0 --n-periods 20000 --min-period 0.3 --max-period 30.0
```

### Step 3: Analyze results

```bash
python3.11 python/analyze_candidates.py --input candidates.json --lightcurves data/lightcurves/ --crossmatch --top-n 30
```

### Step 4: Validate candidates (CRITICAL)

```bash
# Rust validation (fast, parallel via Rayon):
./target/release/hunt validate -i candidates.json -l data/lightcurves -o results/

# Python validation (with ExoFOP cross-reference):
python3.11 python/validate_candidates.py
```

Both run the same 5 false positive tests. The Rust version runs tests in parallel and outputs JSON. The Python version can additionally fetch ExoFOP catalog data for period agreement.

### Step 5: Cross-match (optional)

```bash
./target/release/hunt crossmatch -i candidates.json -c data/lightcurves/confirmed_exoplanets.csv -o results/crossmatch_results.csv
```

## Validation Pipeline (Rust: `validate.rs` / Python: `validate_candidates.py`)

This is the scientific rigor layer. Based on standard methods from:

- **Kovacs, Zucker & Mazeh (2002)** — Original BLS algorithm (A&A, 391, 369-377)
- **NASA SPOC Data Validation** — TESS pipeline false positive diagnostics
- **ExoMiner, LEO-Vetter** — Automated vetting systems used by professionals

### Five False Positive Tests:

1. **Odd/Even Transit Depth Test** — Compares depths of odd vs even transits. Eclipsing binaries show alternating depths (primary/secondary eclipse). A true planet produces identical depths. Pass if depths match within 3-sigma.

2. **Secondary Eclipse Test** — Searches for a brightness dip at orbital phase 0.5. Self-luminous companions (stars, brown dwarfs) produce secondary eclipses. True planets typically do not. Pass if no significant dip at phase 0.5.

3. **Transit Shape Test** — Measures V-shape vs U-shape morphology. Planets produce flat-bottomed (U-shaped) transits due to limb darkening. Grazing binaries produce V-shaped dips. Computed as ratio of ingress slope to flat bottom duration.

4. **Period Agreement Test** — Compares our BLS-detected period to the TESS pipeline period from ExoFOP. Independent period detection at the same value strongly validates the signal. Pass if periods match within 1%.

5. **Radius Ratio Sanity Check** — Computes Rp/Rs from transit depth. Values > 0.3 (planet radius > 30% of star) are physically implausible for planets and indicate an eclipsing binary or blended source.

### Scoring System:

- Each test contributes to a 0-100 **Planet Likelihood Score**
- **High confidence (>=70):** Passes 4+ tests, period matches, small Rp/Rs
- **Medium (50-69):** Passes most tests but some concerns
- **Low (<50):** Multiple failures, likely false positive

### Current Results (200 TOI run):

- 17 high-confidence candidates (score >= 70)
- 114 medium-confidence (50-69)
- 66 low-confidence (<50)
- Most promising: TOI 133.01 (Rp=1.9 R_Earth), TOI 210.01 (Rp=2.2 R_Earth)
- Pipeline validated: TOI 125.04 (already confirmed planet) correctly scored high

## Scientific Methodology & References

### BLS Algorithm

- **Paper:** Kovacs, Zucker & Mazeh (2002), "A box-fitting algorithm in the search for periodic transits", A&A, 391, 369-377
- **ADS:** https://ui.adsabs.harvard.edu/abs/2002A%26A...391..369K
- **Method:** Phase-fold light curve at trial periods, fit box-shaped transit model, compute BLS power statistic
- **Our implementation:** 15,000 log-spaced trial periods (0.5–20 days), 200 phase bins, 8 box widths, SNR >= 6.0 threshold
- **Detection threshold matches literature:** SNR > 6 is standard (Kovacs+2002), NASA SPOC uses 7.1-sigma MES

### NASA SPOC Pipeline Comparison

Our validation tests mirror the NASA SPOC Data Validation (DV) component:

- SPOC runs odd-even depth test — we do this
- SPOC runs centroid offset analysis — we cannot do this (requires pixel-level data)
- SPOC runs model-shift test — we approximate via secondary eclipse search
- SPOC runs background flux analysis — not available from lightkurve downloads
- **Reference:** Twicken et al. (2018, 2019), AAS presentations on TESS SPOC DV

### False Positive Rates in Transit Surveys

- Professional automated systems achieve ~91% true positive, ~97% false positive rejection (LEO-Vetter)
- Our pipeline with 5 tests reduces 197 raw detections to 17 high-confidence (~91% rejection rate)
- Rp/Rs > 1.0 candidates are expected in blind searches and are correctly flagged

### Citizen Science Precedent

- **Planet Hunters TESS:** 284 candidates, 15 confirmed planets (as of 2023)
- **Notable discoveries:** TOI-813 b (Saturn-sized), TOI-4633 c "Percival" (habitable zone, binary system), TOI-6692
- **Our approach is valid** — citizen scientists using BLS detection have published in peer-reviewed journals
- **Reference:** https://heasarc.gsfc.nasa.gov/docs/tess/citizenscience.html

### Submission Path

To submit candidates to NASA:

1. Package candidate with: phase-folded light curve, BLS periodogram, validation test results, stellar parameters
2. Submit as **community TOI (cTOI)** to ExoFOP-TESS: https://exofop.ipac.caltech.edu
3. TESS TOI Team reviews and assigns TOI number if criteria met
4. Community follow-up observations determine final disposition

## Key Commands for Claude Code

**"Hunt for planets" / "Run the pipeline"** → Execute Steps 1-4 in sequence

**"Analyze a specific star"** → Download individually:

```bash
python3.11 -c "
import lightkurve as lk
lc = lk.search_lightcurve('TIC 261136679', mission='TESS').download()
lc = lc.remove_nans().remove_outliers().normalize()
lc.to_csv('data/lightcurves/custom_target.csv')
"
./target/release/hunt -i data/lightcurves -o candidates.json
```

**"Make it viral" / "Shareable content"** →

- Generate best phase-folded plots (top 30 candidates)
- Write compelling REPORT.md with honest, backed-by-data claims
- Update docs/index.html for GitHub Pages
- NEVER fake claims — all results must be reproducible and properly caveated

**"Validate results"** → Run validate_candidates.py and inspect VALIDATION_REPORT.md

## Critical Rules

1. **NEVER claim planet discovery** — Our candidates are "transit detections" requiring independent confirmation
2. **All claims must be backed by peer-reviewed methodology** — cite Kovacs+2002 for BLS, standard DV tests for validation
3. **Rp/Rs > 1.0 means likely binary** — flag these honestly, don't hide them
4. **Use python3.11** — not python3 (which is 3.9 on this system)
5. **Reproducibility matters** — BLS results are deterministic; running 4x gives same top candidates within ~1% SNR
6. **Cross-reference ExoFOP** before making any claims about a specific target

## Dependencies

- **Rust** 1.75+ with: rayon, serde, serde_json, csv, clap, indicatif, anyhow, ordered-float
- **Python** 3.11+ with: lightkurve, astroquery, pandas, numpy, matplotlib, tqdm, scipy
