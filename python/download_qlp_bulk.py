#!/usr/bin/env python3.11
"""Phase 2 FFI Approach 2: Download QLP (Quick Look Pipeline) light curves from MAST.

QLP extracts light curves from TESS Full Frame Images for ~160k-millions of targets
per sector. These stars were NOT on the 2-minute SPOC target list and were NOT
transit-searched by SPOC TPS — making them prime targets for new planet discovery.

CRITICAL: QLP light curves are HLSP products. You MUST query with:
    obs_collection="HLSP"  (NOT "TESS")
    provenance_name="QLP"

This was why earlier attempts returned 0 results.

QLP FITS files contain KSPSAP_FLUX (detrended flux) — the best column for transit search.
TESS-SPOC FFI products use PDCSAP_FLUX instead.

Usage:
  python3.11 python/download_qlp_bulk.py --sector 56 --limit 500
  python3.11 python/download_qlp_bulk.py --sector 56 --limit 2000 --fgk-only
  python3.11 python/download_qlp_bulk.py --sector 56 --limit 500 --author TESS-SPOC
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


def query_qlp_targets(sector: int, author: str = "QLP") -> list[dict]:
    """Query MAST for QLP/TESS-SPOC FFI targets in a sector.

    CRITICAL: Use obs_collection="HLSP" for QLP and TESS-SPOC FFI products.
    Using "TESS" returns 0 results for these HLSP pipelines.

    NOTE: The full HLSP sector query can be very slow (100k+ results).
    If it times out, falls back to querying TESS targets in the sector
    and checking QLP availability per-star via lightkurve.

    Returns list of dicts with tic_id.
    """
    from astroquery.mast import Observations

    # First try the direct HLSP query (fastest if it works)
    print(f"  Querying MAST HLSP for sector {sector} ({author})...")
    print(f"  (This may take a while for QLP — 100k+ targets per sector)")
    try:
        obs = Observations.query_criteria(
            obs_collection="HLSP",
            sequence_number=sector,
            provenance_name=author,
            dataproduct_type="timeseries",
        )

        targets = []
        seen = set()
        for row in obs:
            name = str(row["target_name"]).strip()
            tic_id = None
            if name.isdigit():
                tic_id = int(name)
            elif "TIC" in name.upper():
                import re
                m = re.search(r'(\d{6,})', name)
                if m:
                    tic_id = int(m.group(1))

            if tic_id and tic_id not in seen:
                seen.add(tic_id)
                targets.append({"tic_id": tic_id})

        print(f"  Found {len(targets)} unique {author} targets in sector {sector}")
        return targets

    except Exception as e:
        print(f"  HLSP bulk query failed or timed out: {e}")
        print(f"  Falling back to TESS sector query + per-star QLP check...")
        return _fallback_query_qlp(sector, author)


def _fallback_query_qlp(sector: int, author: str) -> list[dict]:
    """Fallback: get all TESS targets in sector, then check QLP availability.

    This is slower per-star but avoids the massive HLSP sector query.
    Gets stars from the TESS obs_collection (which includes SPOC targets),
    then we check each via lightkurve for QLP availability.
    """
    from astroquery.mast import Observations

    print(f"  Querying TESS targets in sector {sector}...")
    obs = Observations.query_criteria(
        obs_collection="TESS",
        sequence_number=sector,
        dataproduct_type="timeseries",
    )

    tic_ids = set()
    for row in obs:
        name = str(row["target_name"]).strip()
        if name.isdigit():
            tic_ids.add(int(name))

    print(f"  Found {len(tic_ids)} TESS targets — these will be checked for {author} data")
    return [{"tic_id": t} for t in sorted(tic_ids)]


def download_qlp_lightcurve(
    tic_id: int, sector: int, author: str, output_dir: Path
) -> bool:
    """Download a single QLP/TESS-SPOC light curve and convert to CSV.

    Uses lightkurve with author="QLP" or "TESS-SPOC" which handles HLSP lookup.
    Falls back to direct FITS download if lightkurve fails.

    Returns True on success, False on failure.
    """
    import lightkurve as lk

    filename = f"TIC_{tic_id}_s{sector:04d}_qlp.csv"
    filepath = output_dir / filename

    if filepath.exists():
        return True  # Already have it

    try:
        # lightkurve handles HLSP products correctly with author parameter
        search = lk.search_lightcurve(
            f"TIC {tic_id}",
            mission="TESS",
            author=author,
            sector=sector,
        )

        if len(search) == 0:
            return False

        lc = search[0].download(quality_bitmask="hardest")
        if lc is None:
            return False

        lc = lc.remove_nans().remove_outliers(sigma=5).normalize()

        if len(lc.time.value) < 100:
            return False

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
        return True

    except Exception:
        return False


def download_qlp_lightcurve_fits(
    tic_id: int, sector: int, author: str, output_dir: Path
) -> bool:
    """Fallback: download QLP FITS directly from MAST and extract flux.

    QLP uses KSPSAP_FLUX column. TESS-SPOC uses PDCSAP_FLUX.
    """
    from astropy.io import fits
    from astroquery.mast import Observations

    filename = f"TIC_{tic_id}_s{sector:04d}_qlp.csv"
    filepath = output_dir / filename

    if filepath.exists():
        return True

    try:
        obs = Observations.query_criteria(
            obs_collection="HLSP",
            sequence_number=sector,
            provenance_name=author,
            target_name=str(tic_id),
            dataproduct_type="timeseries",
        )

        if len(obs) == 0:
            return False

        products = Observations.get_product_list(obs[0])
        lc_products = products[products["productSubGroupDescription"] == "LLC"]
        if len(lc_products) == 0:
            # Try any FITS product
            lc_products = products[
                [".fits" in str(p).lower() for p in products["productFilename"]]
            ]
        if len(lc_products) == 0:
            return False

        manifest = Observations.download_products(
            lc_products[0], download_dir=str(output_dir / "_fits_cache")
        )
        fits_path = manifest["Local Path"][0]

        with fits.open(fits_path) as hdul:
            data = hdul[1].data
            time = data["TIME"]

            # QLP uses KSPSAP_FLUX, TESS-SPOC uses PDCSAP_FLUX
            flux_col = "KSPSAP_FLUX" if author == "QLP" else "PDCSAP_FLUX"
            if flux_col not in data.columns.names:
                # Fallback to SAP_FLUX
                flux_col = "SAP_FLUX"
            flux = data[flux_col]

            err_col = flux_col.replace("FLUX", "FLUX_ERR")
            if err_col in data.columns.names:
                flux_err = data[err_col]
            else:
                flux_err = np.full_like(flux, 0.001)

        # Clean and normalize
        good = np.isfinite(time) & np.isfinite(flux) & (flux > 0)
        time, flux, flux_err = time[good], flux[good], flux_err[good]

        if len(time) < 100:
            return False

        median_flux = np.median(flux)
        norm_flux = flux / median_flux
        norm_err = flux_err / median_flux

        # Remove 5-sigma outliers
        med = np.median(norm_flux)
        mad = np.median(np.abs(norm_flux - med)) * 1.4826
        keep = np.abs(norm_flux - med) < 5 * mad

        df = pd.DataFrame({
            "time": time[keep],
            "flux": norm_flux[keep],
            "flux_err": norm_err[keep],
        })
        df.to_csv(filepath, index=False)
        return True

    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 FFI: Download QLP/TESS-SPOC bulk light curves from MAST"
    )
    parser.add_argument("--sector", type=int, required=True, help="TESS sector number")
    parser.add_argument("--limit", type=int, default=500, help="Max stars to download")
    parser.add_argument(
        "--author", default="QLP", choices=["QLP", "TESS-SPOC"],
        help="HLSP pipeline (QLP=detrended FFI, TESS-SPOC=PDC FFI)"
    )
    parser.add_argument("--output", default=None, help="Output directory")
    parser.add_argument("--fgk-only", action="store_true", help="Filter for FGK dwarfs")
    parser.add_argument(
        "--use-fits-fallback", action="store_true",
        help="Use direct FITS download instead of lightkurve (slower but more reliable)"
    )
    args = parser.parse_args()

    output_dir = Path(args.output) if args.output else Path(f"data/phase2/qlp_sector_{args.sector}")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"PHASE 2 FFI: QLP/TESS-SPOC Bulk Download — Sector {args.sector}")
    print(f"  Author: {args.author} (obs_collection=HLSP)")
    print(f"  Limit: {args.limit}")
    print(f"  FGK filter: {args.fgk_only}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    # Step 1: Query all QLP targets in sector
    print(f"\n[1/4] Querying MAST for {args.author} targets...")
    targets = query_qlp_targets(args.sector, args.author)

    if not targets:
        print(f"ERROR: No {args.author} targets found for sector {args.sector}.")
        print("  Note: QLP coverage varies by sector. Try a different sector.")
        sys.exit(1)

    tic_ids = [t["tic_id"] for t in targets]
    print(f"  Total {args.author} targets: {len(tic_ids)}")

    # Step 2: Remove known TOIs
    print(f"\n[2/4] Removing known TOIs and cTOIs...")
    sys.path.insert(0, str(Path(__file__).parent))
    from download_sector_bulk import fetch_toi_tic_ids
    known_tics = fetch_toi_tic_ids()
    tic_ids_filtered = [t for t in tic_ids if t not in known_tics]
    print(f"  Removed {len(tic_ids) - len(tic_ids_filtered)} known TOIs/cTOIs")
    print(f"  Remaining: {len(tic_ids_filtered)}")

    # Also remove SPOC 2-minute targets (already transit-searched)
    from download_ffi_tesscut import get_spoc_target_ids
    spoc_tics = get_spoc_target_ids(args.sector)
    before = len(tic_ids_filtered)
    tic_ids_filtered = [t for t in tic_ids_filtered if t not in spoc_tics]
    print(f"  Removed {before - len(tic_ids_filtered)} SPOC 2-min targets")
    print(f"  Non-SPOC, non-TOI targets: {len(tic_ids_filtered)}")

    # Step 3: Optional FGK filter
    if args.fgk_only:
        print(f"\n[3/4] Filtering for FGK dwarfs...")
        from download_sector_bulk import filter_fgk_dwarfs
        tic_ids_filtered = filter_fgk_dwarfs(tic_ids_filtered)
    else:
        print(f"\n[3/4] Skipping FGK filter")

    tic_ids_filtered = tic_ids_filtered[:args.limit]
    print(f"\n  Final target count: {len(tic_ids_filtered)}")

    # Save manifest
    manifest = output_dir / "targets_qlp.json"
    with open(manifest, "w") as f:
        json.dump({
            "sector": args.sector,
            "author": args.author,
            "method": "QLP_HLSP",
            "obs_collection": "HLSP",
            "fgk_only": args.fgk_only,
            "n_targets": len(tic_ids_filtered),
            "tic_ids": tic_ids_filtered,
        }, f, indent=2)

    # Step 4: Download light curves
    print(f"\n[4/4] Downloading {len(tic_ids_filtered)} {args.author} light curves...")
    downloaded = 0
    skipped = 0
    failed = 0

    download_fn = download_qlp_lightcurve_fits if args.use_fits_fallback else download_qlp_lightcurve

    for tic_id in tqdm(tic_ids_filtered, desc=f"  {args.author} download"):
        filepath = output_dir / f"TIC_{tic_id}_s{args.sector:04d}_qlp.csv"
        if filepath.exists():
            skipped += 1
            continue

        success = download_fn(tic_id, args.sector, args.author, output_dir)
        if success:
            downloaded += 1
        else:
            failed += 1

    print(f"\n{'=' * 60}")
    print(f"  QLP DOWNLOAD COMPLETE")
    print(f"  Downloaded:  {downloaded}")
    print(f"  Already had: {skipped}")
    print(f"  Failed:      {failed}")
    print(f"  Output:      {output_dir}")
    print(f"{'=' * 60}")

    if downloaded + skipped > 0:
        print(f"\nNext: Run BLS on QLP light curves:")
        print(f"  ./target/release/hunt search -i {output_dir} -o results/phase2/candidates_qlp_s{args.sector}.json --snr-threshold 5.5 --n-periods 15000")


if __name__ == "__main__":
    main()
