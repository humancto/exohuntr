#!/usr/bin/env python3.11
"""Phase 2 FFI Approach 3: Multi-sector stacking for sub-threshold planet detection.

SPOC TPS has a detection threshold of ~7.1 sigma MES (Multiple Event Statistic).
Our BLS uses SNR >= 6.0. By stacking data from multiple TESS sectors for the same
star, we increase the number of transits observed and push the effective SNR higher.

This approach finds:
  - Long-period planets (P > 15d) that only transit 1-2 times per sector
  - Shallow transits on faint stars that are individually below threshold
  - Planets around stars in the TESS Continuous Viewing Zone (CVZ) with 13+ sectors

Strategy:
  1. Pick a set of TIC IDs (from previous phase2 runs or a target list)
  2. For each star, find ALL available TESS sectors
  3. Download and stitch light curves from all sectors
  4. Save the combined time series for BLS analysis
  5. The longer baseline dramatically improves period determination and SNR

Usage:
  python3.11 python/stack_multisector.py --tic-ids 97168477 14179859 --author SPOC
  python3.11 python/stack_multisector.py --from-json results/phase2/candidates_s56.json --min-snr 5.0
  python3.11 python/stack_multisector.py --sector 56 --top-n 50 --author QLP
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore")


def find_all_sectors(tic_id: int, author: str = "SPOC") -> list[int]:
    """Find all TESS sectors that observed a given TIC ID.

    Queries MAST for all observations of this target across all sectors.
    For HLSP products (QLP, TESS-SPOC), uses obs_collection="HLSP".
    """
    from astroquery.mast import Observations

    obs_collection = "HLSP" if author in ("QLP", "TESS-SPOC") else "TESS"

    try:
        obs = Observations.query_criteria(
            obs_collection=obs_collection,
            provenance_name=author,
            target_name=str(tic_id),
            dataproduct_type="timeseries",
        )

        sectors = set()
        for row in obs:
            seq = row.get("sequence_number")
            if seq is not None:
                try:
                    sectors.add(int(seq))
                except (ValueError, TypeError):
                    pass

        return sorted(sectors)
    except Exception:
        return []


def download_and_stitch(
    tic_id: int, sectors: list[int], author: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]] | None:
    """Download light curves from multiple sectors and stitch them together.

    Each sector's light curve is independently normalized to median=1.0 before
    stitching. This removes sector-to-sector flux offsets from different
    apertures, backgrounds, and CCD positions.

    Returns (time, flux, flux_err, sectors_used) or None on failure.
    """
    import lightkurve as lk

    all_time = []
    all_flux = []
    all_err = []
    sectors_used = []

    for sector in sectors:
        try:
            search = lk.search_lightcurve(
                f"TIC {tic_id}",
                mission="TESS",
                author=author,
                sector=sector,
            )
            if len(search) == 0:
                continue

            lc = search[0].download(quality_bitmask="hardest")
            if lc is None:
                continue

            lc = lc.remove_nans().remove_outliers(sigma=5).normalize()

            if len(lc.time.value) < 50:
                continue

            flux_err = (
                lc.flux_err.value
                if lc.flux_err is not None
                else np.full(len(lc.time.value), 0.001)
            )

            all_time.append(lc.time.value)
            all_flux.append(lc.flux.value)
            all_err.append(flux_err)
            sectors_used.append(sector)

        except Exception:
            continue

    if not all_time:
        return None

    time = np.concatenate(all_time)
    flux = np.concatenate(all_flux)
    flux_err = np.concatenate(all_err)

    # Sort by time
    order = np.argsort(time)
    time = time[order]
    flux = flux[order]
    flux_err = flux_err[order]

    if len(time) < 200:
        return None

    return time, flux, flux_err, sectors_used


def get_candidates_from_json(json_path: str, min_snr: float = 0.0) -> list[int]:
    """Extract TIC IDs from a BLS candidates JSON file.

    Optionally filter by minimum SNR to focus on marginal detections
    that might become significant with more data.
    """
    with open(json_path) as f:
        candidates = json.load(f)

    tic_ids = []
    for c in candidates:
        snr = c.get("snr", 0)
        if snr >= min_snr:
            filename = c.get("filename", "")
            # Extract TIC ID from filename like "TIC_12345678_s0056.csv"
            import re
            m = re.search(r'TIC_(\d+)', filename)
            if m:
                tic_ids.append(int(m.group(1)))

    return tic_ids


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 FFI: Multi-sector stacking for sub-threshold planet detection"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--tic-ids", nargs="+", type=int, help="Specific TIC IDs to stack")
    group.add_argument("--from-json", type=str, help="BLS candidates JSON file")
    group.add_argument("--sector", type=int, help="Re-stack all candidates from a sector run")

    parser.add_argument("--min-snr", type=float, default=0.0,
                        help="Min SNR filter when using --from-json (default: all)")
    parser.add_argument("--top-n", type=int, default=50,
                        help="Max targets to stack (default: 50)")
    parser.add_argument(
        "--author", default="SPOC", choices=["SPOC", "QLP", "TESS-SPOC"],
        help="Light curve pipeline to use for download"
    )
    parser.add_argument("--output", default=None, help="Output directory")
    parser.add_argument("--min-sectors", type=int, default=2,
                        help="Minimum sectors required (default: 2)")
    args = parser.parse_args()

    # Determine TIC IDs
    if args.tic_ids:
        tic_ids = args.tic_ids
    elif args.from_json:
        tic_ids = get_candidates_from_json(args.from_json, args.min_snr)
    else:
        # Load from sector run
        json_path = f"results/phase2/candidates_s{args.sector}.json"
        if not Path(json_path).exists():
            print(f"ERROR: {json_path} not found. Run phase2-hunt first.")
            sys.exit(1)
        tic_ids = get_candidates_from_json(json_path, args.min_snr)

    tic_ids = tic_ids[:args.top_n]

    output_dir = Path(args.output) if args.output else Path("data/phase2/multisector")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"PHASE 2: Multi-Sector Stacking")
    print(f"  Targets: {len(tic_ids)}")
    print(f"  Author: {args.author}")
    print(f"  Min sectors: {args.min_sectors}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    stacked = 0
    skipped_few_sectors = 0
    failed = 0
    results_log = []

    for tic_id in tqdm(tic_ids, desc="  Multi-sector stacking"):
        filename = f"TIC_{tic_id}_multisector.csv"
        filepath = output_dir / filename

        if filepath.exists():
            stacked += 1
            continue

        # Find all available sectors
        sectors = find_all_sectors(tic_id, args.author)

        if len(sectors) < args.min_sectors:
            skipped_few_sectors += 1
            results_log.append({
                "tic_id": tic_id,
                "status": "skipped",
                "reason": f"only {len(sectors)} sectors (need {args.min_sectors})",
                "sectors_available": sectors,
            })
            continue

        # Download and stitch
        result = download_and_stitch(tic_id, sectors, args.author)

        if result is None:
            failed += 1
            results_log.append({
                "tic_id": tic_id,
                "status": "failed",
                "sectors_available": sectors,
            })
            continue

        time, flux, flux_err, sectors_used = result
        baseline_days = time[-1] - time[0]

        df = pd.DataFrame({
            "time": time,
            "flux": flux,
            "flux_err": flux_err,
        })
        df.to_csv(filepath, index=False)

        stacked += 1
        results_log.append({
            "tic_id": tic_id,
            "status": "stacked",
            "sectors_used": sectors_used,
            "n_sectors": len(sectors_used),
            "n_points": len(time),
            "baseline_days": round(baseline_days, 1),
        })

        tqdm.write(
            f"    TIC {tic_id}: {len(sectors_used)} sectors, "
            f"{len(time)} pts, {baseline_days:.0f} days"
        )

    # Save stacking log
    log_path = output_dir / "stacking_log.json"
    with open(log_path, "w") as f:
        json.dump(results_log, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"  MULTI-SECTOR STACKING COMPLETE")
    print(f"  Stacked:            {stacked}")
    print(f"  Too few sectors:    {skipped_few_sectors}")
    print(f"  Failed:             {failed}")
    print(f"  Output:             {output_dir}")
    print(f"  Log:                {log_path}")
    print(f"{'=' * 60}")

    # Summary of best targets
    stacked_results = [r for r in results_log if r.get("status") == "stacked"]
    if stacked_results:
        stacked_results.sort(key=lambda r: r.get("n_sectors", 0), reverse=True)
        print(f"\nTop stacked targets (by sector count):")
        for r in stacked_results[:10]:
            print(
                f"  TIC {r['tic_id']}: {r['n_sectors']} sectors, "
                f"{r['n_points']} pts, {r['baseline_days']} days"
            )

    if stacked > 0:
        print(f"\nNext: Run BLS on stacked light curves:")
        print(f"  ./target/release/hunt search -i {output_dir} -o results/phase2/candidates_multisector.json --snr-threshold 5.5 --n-periods 20000 --max-period 40.0")
        print(f"\n  NOTE: Use --max-period 40.0 to search for long-period planets")
        print(f"  that single-sector BLS would miss!")


if __name__ == "__main__":
    main()
