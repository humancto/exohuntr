#!/usr/bin/env python3.11
"""Download light curves from less-studied TESS sectors (80-96) for new detections.

Strategy: Query ExoFOP TOI list, filter for TOIs observed in sectors 80-96,
and download those that are still unconfirmed (PC disposition).
"""
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    import os
    import sys
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from tqdm import tqdm

    OUTPUT_DIR = Path('data/lightcurves_new_sectors')
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    SECTOR_MIN = 80
    SECTOR_MAX = 96
    LIMIT = 200

    print("=" * 60)
    print(f"Downloading TOIs from TESS Sectors {SECTOR_MIN}-{SECTOR_MAX}")
    print("=" * 60)

    # Step 1: Get TOI catalog from ExoFOP
    print("\n[1/3] Fetching TOI catalog from ExoFOP...", flush=True)
    try:
        toi_url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        tois = pd.read_csv(toi_url, comment="#")
        print(f"  Total TOIs in catalog: {len(tois)}", flush=True)
    except Exception as e:
        print(f"  ERROR: Could not fetch TOI list: {e}", flush=True)
        sys.exit(1)

    # Step 2: Filter for unconfirmed candidates in target sectors
    print(f"\n[2/3] Filtering for unconfirmed TOIs in sectors {SECTOR_MIN}-{SECTOR_MAX}...", flush=True)

    # Keep only Planet Candidates
    candidates = tois[tois["TFOPWG Disposition"].isin(["PC", ""])]
    print(f"  Unconfirmed candidates: {len(candidates)}", flush=True)

    # Filter by sector — the "Sectors" column contains comma-separated sector numbers
    def in_target_sectors(sectors_str):
        try:
            sectors = [int(s.strip()) for s in str(sectors_str).split(',')]
            return any(SECTOR_MIN <= s <= SECTOR_MAX for s in sectors)
        except (ValueError, AttributeError):
            return False

    if 'Sectors' in candidates.columns:
        target_candidates = candidates[candidates['Sectors'].apply(in_target_sectors)]
    elif 'Sector' in candidates.columns:
        target_candidates = candidates[candidates['Sector'].apply(in_target_sectors)]
    else:
        # Try to find the right column
        print(f"  Available columns: {list(candidates.columns)}", flush=True)
        print("  WARNING: No 'Sectors' column found. Downloading general unconfirmed TOIs.", flush=True)
        # Fall back: take TOIs with high TOI numbers (newer, likely from later sectors)
        candidates_sorted = candidates.sort_values('TOI', ascending=False)
        target_candidates = candidates_sorted.head(LIMIT)

    target_candidates = target_candidates.head(LIMIT)
    print(f"  TOIs in sectors {SECTOR_MIN}-{SECTOR_MAX}: {len(target_candidates)}", flush=True)

    if len(target_candidates) == 0:
        print("  No TOIs found in target sectors. Trying newest unconfirmed TOIs instead...", flush=True)
        candidates_sorted = candidates.sort_values('TOI', ascending=False)
        target_candidates = candidates_sorted.head(LIMIT)
        print(f"  Using {len(target_candidates)} newest unconfirmed TOIs", flush=True)

    # Step 3: Download light curves
    print(f"\n[3/3] Downloading {len(target_candidates)} light curves...", flush=True)
    import lightkurve as lk

    downloaded = 0
    failed = 0
    already_have = 0

    for _, row in tqdm(target_candidates.iterrows(), total=len(target_candidates), desc="  Downloading"):
        try:
            tic_id = f"TIC {int(row['TIC ID'])}"
            toi_num = str(row.get('TOI', 'unknown'))

            # Check if we already have this
            filename = f"TOI_{toi_num}_{tic_id.replace(' ', '_')}.csv"
            filepath = OUTPUT_DIR / filename

            # Also check the original data dir
            orig_path = Path('data/lightcurves') / filename
            if filepath.exists() or orig_path.exists():
                already_have += 1
                continue

            search = lk.search_lightcurve(tic_id, mission='TESS', author='SPOC')
            if len(search) == 0:
                failed += 1
                continue

            lc = search[0].download(quality_bitmask='hardest')
            if lc is None:
                failed += 1
                continue

            lc = lc.remove_nans().remove_outliers(sigma=5).normalize()
            if len(lc.time.value) < 100:
                failed += 1
                continue

            flux_err = lc.flux_err.value if lc.flux_err is not None else np.full(len(lc.time.value), 0.001)
            df = pd.DataFrame({
                'time': lc.time.value,
                'flux': lc.flux.value,
                'flux_err': flux_err,
            })
            df.to_csv(filepath, index=False)
            downloaded += 1

        except Exception as e:
            failed += 1
            continue

    print(f"\n  Downloaded: {downloaded}")
    print(f"  Already had: {already_have}")
    print(f"  Failed: {failed}")
    print(f"  Output: {OUTPUT_DIR}/")

    if downloaded > 0:
        print(f"\n  Next: Run BLS on these:")
        print(f"  cargo build --release")
        print(f"  ./target/release/hunt search -i {OUTPUT_DIR} -o candidates_new_sectors.json --snr-threshold 6.0")
        print(f"  ./target/release/hunt validate -i candidates_new_sectors.json -l {OUTPUT_DIR} -o results/")

    print("\nDone!", flush=True)
