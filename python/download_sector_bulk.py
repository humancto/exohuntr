#!/usr/bin/env python3.11
"""Phase 2: Download ALL light curves from a TESS sector for new planet discovery.

Unlike Phase 1 (which only downloads known TOIs), this script downloads light
curves for ALL observed stars in a sector, then filters out known TOIs so we
can run BLS on unstudied stars. Any transit detection on a non-TOI star is
potentially a NEW planet candidate.

Strategy:
  1. Query MAST for all light curves in a sector (SPOC or QLP)
  2. Download the ExoFOP TOI list and community TOI (cTOI) list
  3. Filter out stars with existing TOI/cTOI designations
  4. Optionally filter for FGK dwarfs (Tmag 8-13) using TIC catalog
  5. Download light curves in parallel batches with resume support
  6. Save as CSV (time, flux, flux_err) — same format as Phase 1

Usage:
  python3.11 python/download_sector_bulk.py --sector 56 --limit 1000
  python3.11 python/download_sector_bulk.py --sector 70 --limit 5000 --author QLP
  python3.11 python/download_sector_bulk.py --sector 56 --limit 2000 --fgk-only
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


def fetch_toi_tic_ids() -> set[int]:
    """Download TOI + cTOI lists from ExoFOP and return all known TIC IDs."""
    known_tics = set()

    # TOI list
    try:
        toi_url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        tois = pd.read_csv(toi_url, comment="#")
        known_tics.update(tois["TIC ID"].dropna().astype(int).unique())
        print(f"  TOI catalog: {len(tois)} entries, {len(known_tics)} unique TIC IDs")
    except Exception as e:
        print(f"  WARNING: Could not fetch TOI list: {e}")

    # Community TOI list
    try:
        ctoi_url = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
        ctois = pd.read_csv(ctoi_url, comment="#")
        ctoi_tics = set(ctois["TIC ID"].dropna().astype(int).unique())
        known_tics.update(ctoi_tics)
        print(f"  cTOI catalog: {len(ctois)} entries, {len(ctoi_tics)} unique TIC IDs")
    except Exception as e:
        print(f"  WARNING: Could not fetch cTOI list: {e}")

    return known_tics


def query_sector_targets(sector: int, author: str = "SPOC") -> list[dict]:
    """Query MAST for all light curve targets in a sector.

    Returns list of dicts with 'target_name' and 'obsid' fields.
    """
    from astroquery.mast import Observations

    print(f"  Querying MAST for sector {sector} ({author})...")
    obs = Observations.query_criteria(
        obs_collection="TESS",
        sequence_number=sector,
        provenance_name=author,
        dataproduct_type="timeseries",
    )

    targets = []
    seen = set()
    for row in obs:
        name = str(row["target_name"]).strip()
        if name.isdigit() and name not in seen:
            seen.add(name)
            targets.append({
                "tic_id": int(name),
                "obsid": str(row["obsid"]),
            })

    print(f"  Found {len(targets)} unique targets in sector {sector}")
    return targets


def filter_fgk_dwarfs(tic_ids: list[int], batch_size: int = 500) -> list[int]:
    """Filter TIC IDs to FGK dwarf stars (best for transit detection).

    Criteria: Tmag 8-13, logg > 4.0, Teff 3500-7000 K
    """
    from astroquery.mast import Catalogs

    print(f"  Filtering {len(tic_ids)} stars for FGK dwarfs...")
    good = []

    for i in tqdm(range(0, len(tic_ids), batch_size), desc="  TIC query"):
        batch = tic_ids[i:i + batch_size]
        try:
            result = Catalogs.query_criteria(
                catalog="TIC",
                ID=batch,
            )
            for row in result:
                tmag = row.get("Tmag")
                logg = row.get("logg")
                teff = row.get("Teff")

                if (tmag is not None and 8 < float(tmag) < 13
                        and logg is not None and float(logg) > 4.0
                        and teff is not None and 3500 < float(teff) < 7000):
                    good.append(int(row["ID"]))
        except Exception as e:
            print(f"  WARNING: TIC batch query failed: {e}")
            # On failure, include the batch unfiltered
            good.extend(batch)

    print(f"  FGK dwarfs: {len(good)} / {len(tic_ids)}")
    return good


def download_lightcurves(
    tic_ids: list[int],
    sector: int,
    output_dir: Path,
    author: str = "SPOC",
) -> tuple[int, int, int]:
    """Download light curves for a list of TIC IDs.

    Returns (downloaded, skipped_existing, failed).
    """
    import lightkurve as lk

    downloaded = 0
    skipped = 0
    failed = 0

    for tic_id in tqdm(tic_ids, desc="  Downloading"):
        filename = f"TIC_{tic_id}_s{sector:04d}.csv"
        filepath = output_dir / filename

        if filepath.exists():
            skipped += 1
            continue

        try:
            search = lk.search_lightcurve(
                f"TIC {tic_id}",
                mission="TESS",
                author=author,
                sector=sector,
            )
            if len(search) == 0:
                failed += 1
                continue

            lc = search[0].download(quality_bitmask="hardest")
            if lc is None:
                failed += 1
                continue

            lc = lc.remove_nans().remove_outliers(sigma=5).normalize()

            if len(lc.time.value) < 100:
                failed += 1
                continue

            flux_err = (
                lc.flux_err.value
                if lc.flux_err is not None
                else np.full(len(lc.time.value), 0.001)
            )

            df = pd.DataFrame({
                "time": lc.time.value,
                "flux": lc.flux.value,
                "flux_err": flux_err,
            })
            df.to_csv(filepath, index=False)
            downloaded += 1

        except Exception:
            failed += 1
            continue

    return downloaded, skipped, failed


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2: Bulk sector download for new planet discovery"
    )
    parser.add_argument("--sector", type=int, required=True, help="TESS sector number")
    parser.add_argument("--limit", type=int, default=1000, help="Max stars to download")
    parser.add_argument(
        "--author", default="SPOC", choices=["SPOC", "QLP", "TESS-SPOC"],
        help="Light curve pipeline (SPOC=2min, QLP=FFI, TESS-SPOC=FFI)"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output directory (default: data/phase2/sector_N)"
    )
    parser.add_argument(
        "--fgk-only", action="store_true",
        help="Filter for FGK dwarf stars (Tmag 8-13, logg>4, Teff 3500-7000K)"
    )
    parser.add_argument(
        "--include-tois", action="store_true",
        help="Include known TOIs (default: exclude them for new discovery)"
    )
    args = parser.parse_args()

    output_dir = Path(args.output) if args.output else Path(f"data/phase2/sector_{args.sector}")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"PHASE 2: Bulk Sector Download — Sector {args.sector}")
    print(f"  Author: {args.author}")
    print(f"  Limit: {args.limit}")
    print(f"  FGK filter: {args.fgk_only}")
    print(f"  Include TOIs: {args.include_tois}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    # Step 1: Query all targets in sector
    print(f"\n[1/4] Querying MAST for sector {args.sector} targets...")
    targets = query_sector_targets(args.sector, args.author)

    if not targets:
        print("ERROR: No targets found. Check sector number and author.")
        sys.exit(1)

    tic_ids = [t["tic_id"] for t in targets]

    # Step 2: Filter out known TOIs
    if not args.include_tois:
        print("\n[2/4] Filtering out known TOIs and cTOIs...")
        known_tics = fetch_toi_tic_ids()
        before = len(tic_ids)
        tic_ids = [t for t in tic_ids if t not in known_tics]
        print(f"  Removed {before - len(tic_ids)} known TOI/cTOI stars")
        print(f"  Remaining non-TOI targets: {len(tic_ids)}")
    else:
        print("\n[2/4] Skipping TOI filter (--include-tois)")

    # Step 3: Optional FGK filter
    if args.fgk_only:
        print("\n[3/4] Filtering for FGK dwarf stars...")
        tic_ids = filter_fgk_dwarfs(tic_ids)
    else:
        print("\n[3/4] Skipping FGK filter (use --fgk-only to enable)")

    # Apply limit
    tic_ids = tic_ids[:args.limit]
    print(f"\n  Final target count: {len(tic_ids)}")

    # Save target list for reproducibility
    manifest = output_dir / "targets.json"
    with open(manifest, "w") as f:
        json.dump({
            "sector": args.sector,
            "author": args.author,
            "fgk_only": args.fgk_only,
            "include_tois": args.include_tois,
            "n_targets": len(tic_ids),
            "tic_ids": tic_ids,
        }, f, indent=2)
    print(f"  Target list saved: {manifest}")

    # Step 4: Download
    print(f"\n[4/4] Downloading {len(tic_ids)} light curves...")
    downloaded, skipped, failed = download_lightcurves(
        tic_ids, args.sector, output_dir, args.author
    )

    print(f"\n{'=' * 60}")
    print(f"  DOWNLOAD COMPLETE")
    print(f"  Downloaded:  {downloaded}")
    print(f"  Already had: {skipped}")
    print(f"  Failed:      {failed}")
    print(f"  Output:      {output_dir}")
    print(f"{'=' * 60}")

    if downloaded + skipped > 0:
        print(f"\nNext steps:")
        print(f"  1. Run BLS:")
        print(f"     ./target/release/hunt search -i {output_dir} -o results/phase2/candidates_s{args.sector}.json --snr-threshold 6.0 --n-periods 15000")
        print(f"  2. Validate:")
        print(f"     ./target/release/hunt validate -i results/phase2/candidates_s{args.sector}.json -l {output_dir} -o results/phase2/")
        print(f"  3. Flag discoveries:")
        print(f"     python3.11 python/flag_discoveries.py --input results/phase2/candidates_s{args.sector}.json --validation results/phase2/validation_results.json")


if __name__ == "__main__":
    main()
