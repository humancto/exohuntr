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

---

## Run 2: Sector 70, FGK Dwarfs (Partial — 203 stars)

### Parameters
- **Date:** 2026-03-29
- **Sector:** 70 (chosen for lowest TOI coverage: 1.2%)
- **Author:** SPOC
- **FGK filter:** ON (Tmag 8-13, logg>4, Teff 3500-7000K)
- **Stars downloaded:** 203 / 500 target (still downloading)

### Results

| Stage | Count |
|-------|-------|
| Light curves | 203 |
| BLS detections | 176 |
| Validated (score >= 70) | 45 |
| Flagged NEW | 88 |
| Planet-sized NEW (Rp/Rs < 0.3) | 0 |

Smallest Rp/Rs among NEW candidates: TIC 14179859 (Rp/Rs=0.385, P=10.28d, SNR=9.7)

### Conclusion

Same pattern as sector 56: no planet-sized candidates among SPOC non-TOI targets. All detections are eclipsing binaries (Rp/Rs > 0.3). The FGK filter didn't change the outcome — the fundamental issue is that SPOC targets have already been searched by the TESS pipeline.

## Strategic Assessment

### Why SPOC non-TOIs don't yield new planets

The TESS SPOC pipeline runs its own transit search (using the Transiting Planet Search module, TPS) on every 2-minute cadence target. Any planet candidate it finds becomes a Threshold Crossing Event (TCE), goes through Data Validation (DV), and if it passes, becomes a TOI.

**Stars in the SPOC target list that are NOT TOIs have already been searched and cleared by SPOC.** Our BLS is finding the same eclipsing binaries that SPOC found and correctly rejected.

### What would actually find new planets

1. **Full Frame Image (FFI) data** — 200k+ stars per sector that were NOT on the 2-minute target list and were NOT searched by SPOC TPS. These stars only have FFI photometry (10-min or 200-sec cadence). QLP extracts light curves but doesn't run a transit search.

2. **Problem:** QLP/TESS-SPOC FFI light curves are not available as individual timeseries products on MAST through astroquery. They may need to be downloaded as bulk FITS files from the MAST archive directly.

3. **Alternative approach:** Use `eleanor` Python package to extract light curves directly from FFI cutouts for specific TIC IDs. This bypasses the MAST product catalog entirely.

4. **Another alternative:** Download pre-computed QLP light curves from the bulk download portal at https://archive.stsci.edu/hlsp/qlp

### Recommended next steps

- Investigate bulk QLP download from MAST HLSP portal
- Try `eleanor` for on-demand FFI extraction
- Or: accept that laptop-based citizen science may not compete with SPOC on 2-min targets, and focus on:
  - Multi-sector stacking for marginal detections below SPOC's threshold
  - Long-period planets that SPOC misses (P > 15 days, single transit events)
  - Unusual transit shapes that automated pipelines reject
