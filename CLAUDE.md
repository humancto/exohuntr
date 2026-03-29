# Exoplanet Hunter 🔭

## Project Overview
A Rust + Python exoplanet discovery pipeline that downloads TESS/Kepler light curves and runs BLS (Box-fitting Least Squares) transit detection to find planet candidates. The Rust engine handles the compute-heavy BLS search in parallel; Python handles data download, analysis, and visualization.

## Architecture
```
exoplanet-hunter/
├── Cargo.toml              # Rust project config
├── src/main.rs             # Rust BLS engine (the fast part)
├── python/
│   ├── download_lightcurves.py  # Fetch data from MAST/NASA
│   └── analyze_candidates.py    # Analyze results, generate plots
├── data/
│   └── lightcurves/         # Downloaded CSV light curves
├── results/
│   ├── plots/               # Phase-folded light curve PNGs
│   ├── REPORT.md            # Analysis report
│   └── crossmatch_results.csv
├── candidates.json          # BLS output (intermediate)
├── Makefile                 # Build & run automation
└── CLAUDE.md                # You are here
```

## How to Run the Full Pipeline

### Step 1: Download light curves
```bash
cd python
pip install lightkurve astroquery pandas numpy matplotlib tqdm
python download_lightcurves.py --mission tess --sector 56 --limit 300 --catalog
```

For unconfirmed candidates (best chance at discovery):
```bash
python download_lightcurves.py --candidates-only --limit 500 --catalog
```

### Step 2: Build and run the Rust BLS engine
```bash
cargo build --release
./target/release/hunt -i data/lightcurves -o candidates.json --snr-threshold 6.0
```

Aggressive search (more candidates, more false positives):
```bash
./target/release/hunt -i data/lightcurves -o candidates.json \
  --snr-threshold 5.0 --n-periods 20000 --min-period 0.3 --max-period 30.0
```

### Step 3: Analyze results
```bash
python python/analyze_candidates.py --input candidates.json --lightcurves data/lightcurves/ --crossmatch --top-n 30
```

## Key Commands for Claude Code

When asked to "hunt for planets" or "run the pipeline", execute all three steps above in sequence.

When asked to "analyze a specific star", download its light curve individually:
```bash
python -c "
import lightkurve as lk
lc = lk.search_lightcurve('TIC 261136679', mission='TESS').download()
lc = lc.remove_nans().remove_outliers().normalize()
lc.to_csv('data/lightcurves/custom_target.csv')
"
```

When asked to "make it viral" or "make shareable content":
- Generate the best phase-folded plots
- Write a compelling REPORT.md with the top discoveries
- Create a GitHub README with results and pretty plots
- Suggest posting to r/Astronomy, r/exoplanets, or Twitter/X

## Important Scientific Notes
- SNR > 7 is a strong candidate; SNR 5-7 is marginal
- Need at least 2-3 transits for a credible detection
- Always check for: eclipsing binaries (V-shaped dips), centroid shifts, odd/even depth differences
- Cross-reference with ExoFOP (https://exofop.ipac.caltech.edu/tess/) before claiming discovery
- Submit real candidates to AAVSO Exoplanet Database

## Dependencies
- **Rust**: 1.75+ (install via rustup)
- **Python**: 3.10+ with lightkurve, astroquery, pandas, numpy, matplotlib, tqdm
