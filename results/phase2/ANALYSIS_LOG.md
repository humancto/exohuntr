# Phase 2 Analysis Log — Sector 56

## Run Parameters
- **Date:** 2026-03-29
- **Sector:** 56
- **Author:** SPOC (2-min cadence)
- **Stars downloaded:** 883 / 1000 targets (117 failed downloads)
- **TOI filter:** Excluded all known TOIs and cTOIs
- **FGK filter:** Off (all spectral types)

## Pipeline Results

| Stage | Count |
|-------|-------|
| Light curves downloaded | 883 |
| BLS detections (SNR >= 6.0) | 828 |
| Validated (score >= 70) | 274 |
| Validated (score >= 60) | 517 |
| Flagged NEW (not in TOI/cTOI/confirmed) | 434 |
| Planet-sized NEW (Rp/Rs < 0.3) | 1 |

## Candidate Investigation: TIC 97168477

### BLS Detection
- Period: 12.3322 d
- Depth: 76,648 ppm
- SNR: 6.5
- Rp/Rs: 0.2769
- N transits: 2
- Planet score: 80/100

### Stellar Parameters (TIC Catalog)
- Teff: 5383 K (G/K dwarf)
- R*: 0.909 R_sun
- M*: 0.934 M_sun
- logg: 4.49 (main-sequence dwarf)
- Tmag: 11.7
- Distance: 237 pc
- Luminosity class: DWARF

### Derived Planet Parameters (if real)
- Rp = 2.45 R_Jupiter = 27.5 R_Earth (suspiciously large)
- Semi-major axis: 0.102 AU
- T_eq ~ 709 K

### Deep Analysis Verdict: FALSE POSITIVE

**Evidence:**
1. Phase-folded light curve shows NO transit dip. In-transit median flux (1.000051) is actually *higher* than out-of-transit (0.999997). Measured depth = -54 ppm (negative).
2. Multi-sector data (sectors 56 + 83, 759 days, 136,883 points) shows no periodic signal at P=12.33d.
3. Only 2 "transits" in sector 56's 27.9-day baseline — insufficient for reliable detection.
4. Deepest flux point (0.9731) is a 4.4-sigma outlier, not a transit.
5. BLS detected scatter/noise pattern, not a physical transit.
6. 3.0-sigma secondary eclipse hint at phase 0.5 — consistent with systematic noise, not a self-luminous companion.

**Plots:** `TIC_97168477_phase_fold.png`, `TIC_97168477_sector56_detail.png`

## Conclusions

### Sector 56 SPOC Results
- **No genuine new planet candidates found** in 883 non-TOI SPOC targets
- All 434 NEW detections have Rp/Rs > 0.3 (eclipsing binaries) except TIC 97168477
- TIC 97168477 (Rp/Rs=0.277) is confirmed false positive upon deep analysis
- This is expected: SPOC 2-minute targets are a curated list where the TESS pipeline already searched for planets

### Next Steps
1. **Try QLP data** — 10x more stars per sector from Full Frame Images, many never searched
2. **Target recent sectors (70+)** with less community scrutiny
3. **Use --fgk-only filter** to focus on FGK dwarfs (best planet hosts)
4. **Lower SNR threshold cautiously** (5.5 instead of 6.0) for marginal but real signals
