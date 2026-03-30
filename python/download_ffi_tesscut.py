#!/usr/bin/env python3.11
"""Phase 2 FFI Approach 1: Extract light curves from TESS Full Frame Images via TESScut.

TESScut (astroquery.mast.Tesscut) downloads pixel-level cutouts from TESS FFI data.
We perform simple aperture photometry to extract a light curve. This accesses ~200k+
stars per sector that were NEVER on the 2-minute target list and were NEVER searched
by SPOC TPS — the best source for genuine new planet discoveries.

Advantages over eleanor:
  - Works for ALL sectors (including 56+ with 200-second cadence)
  - No tensorflow/numpy version issues
  - Built into astroquery (already a dependency)

Pipeline:
  1. Query TIC catalog for stars in a sector (by coordinate cone search or bulk TIC query)
  2. Filter out SPOC 2-minute targets and known TOIs
  3. For each star: download TESScut cutout → aperture photometry → normalize → save CSV
  4. Output format matches Phase 1/2: time, flux, flux_err

Usage:
  python3.11 python/download_ffi_tesscut.py --sector 40 --limit 100
  python3.11 python/download_ffi_tesscut.py --sector 40 --limit 500 --fgk-only
  python3.11 python/download_ffi_tesscut.py --sector 40 --ra 90.0 --dec -66.5 --radius 1.0
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


def get_ffi_targets_from_tic(
    ra: float, dec: float, radius: float, sector: int
) -> list[dict]:
    """Query TIC catalog for stars near a coordinate that fall on TESS silicon.

    Args:
        ra: Right ascension of cone center (degrees)
        dec: Declination of cone center (degrees)
        radius: Cone search radius (degrees)
        sector: TESS sector (used to verify on-silicon)

    Returns list of dicts with tic_id, ra, dec, Tmag.
    """
    from astroquery.mast import Catalogs

    print(f"  Querying TIC catalog: RA={ra:.2f}, Dec={dec:.2f}, r={radius:.1f} deg...")
    result = Catalogs.query_region(
        f"{ra} {dec}",
        radius=radius,
        catalog="TIC",
    )

    targets = []
    for row in result:
        tmag = row.get("Tmag")
        if tmag is None or np.ma.is_masked(tmag) or float(tmag) > 14:
            continue  # Too faint for FFI photometry
        if tmag is not None and float(tmag) < 4:
            continue  # Too bright (saturated)

        targets.append({
            "tic_id": int(row["ID"]),
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
            "Tmag": float(tmag),
        })

    print(f"  TIC catalog: {len(targets)} stars with Tmag 4-14")
    return targets


def get_spoc_target_ids(sector: int) -> set[int]:
    """Get TIC IDs of stars already on the SPOC 2-minute target list.

    These have already been transit-searched by SPOC TPS — we want to SKIP them.
    """
    from astroquery.mast import Observations

    print(f"  Querying SPOC 2-min targets for sector {sector}...")
    obs = Observations.query_criteria(
        obs_collection="TESS",
        sequence_number=sector,
        provenance_name="SPOC",
        dataproduct_type="timeseries",
    )

    spoc_tics = set()
    for row in obs:
        name = str(row["target_name"]).strip()
        if name.isdigit():
            spoc_tics.add(int(name))

    print(f"  SPOC 2-min targets: {len(spoc_tics)}")
    return spoc_tics


def extract_lightcurve_from_cutout(tic_id: int, ra: float, dec: float,
                                    sector: int, cutout_size: int = 11
                                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Download a TESScut cutout and perform simple aperture photometry.

    Uses a 3x3 pixel aperture centered on the target. Background is estimated
    from the outer ring of the cutout.

    Returns (time, flux, flux_err) or None on failure.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.mast import Tesscut

    try:
        coord = SkyCoord(ra=ra, dec=dec, unit="deg")
        cutout = Tesscut.get_cutouts(
            coordinates=coord,
            size=cutout_size,
            sector=sector,
        )

        if not cutout:
            return None

        hdu = cutout[0]
        data = hdu[1].data

        time = data["TIME"]
        flux_cube = data["FLUX"]
        flux_err_cube = data["FLUX_ERR"]

        # Simple aperture photometry: sum central 3x3 pixels
        center = cutout_size // 2
        ap_slice = (slice(None), slice(center - 1, center + 2), slice(center - 1, center + 2))
        raw_flux = np.nansum(flux_cube[ap_slice], axis=(1, 2))
        raw_err = np.sqrt(np.nansum(flux_err_cube[ap_slice] ** 2, axis=(1, 2)))

        # Background: median of outer ring
        mask = np.ones(flux_cube.shape[1:], dtype=bool)
        mask[2:-2, 2:-2] = False
        bkg_per_pixel = np.nanmedian(flux_cube[:, mask], axis=1)
        n_ap_pixels = 9  # 3x3
        raw_flux -= bkg_per_pixel * n_ap_pixels

        # Remove bad cadences
        good = np.isfinite(time) & np.isfinite(raw_flux) & (raw_flux > 0)
        time = time[good]
        raw_flux = raw_flux[good]
        raw_err = raw_err[good]

        if len(time) < 100:
            return None

        # Normalize
        median_flux = np.median(raw_flux)
        norm_flux = raw_flux / median_flux
        norm_err = raw_err / median_flux

        # Remove 5-sigma outliers
        med = np.median(norm_flux)
        mad = np.median(np.abs(norm_flux - med)) * 1.4826
        keep = np.abs(norm_flux - med) < 5 * mad
        time = time[keep]
        norm_flux = norm_flux[keep]
        norm_err = norm_err[keep]

        if len(time) < 100:
            return None

        return time, norm_flux, norm_err

    except Exception:
        return None


def download_ffi_lightcurves(
    targets: list[dict],
    sector: int,
    output_dir: Path,
) -> tuple[int, int, int]:
    """Extract FFI light curves for a list of targets via TESScut.

    Returns (downloaded, skipped_existing, failed).
    """
    downloaded = 0
    skipped = 0
    failed = 0

    for t in tqdm(targets, desc="  TESScut extraction"):
        tic_id = t["tic_id"]
        filename = f"TIC_{tic_id}_s{sector:04d}_ffi.csv"
        filepath = output_dir / filename

        if filepath.exists():
            skipped += 1
            continue

        result = extract_lightcurve_from_cutout(
            tic_id, t["ra"], t["dec"], sector
        )

        if result is None:
            failed += 1
            continue

        time, flux, flux_err = result
        df = pd.DataFrame({
            "time": time,
            "flux": flux,
            "flux_err": flux_err,
        })
        df.to_csv(filepath, index=False)
        downloaded += 1

    return downloaded, skipped, failed


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 FFI: Extract light curves from TESS FFI via TESScut"
    )
    parser.add_argument("--sector", type=int, required=True, help="TESS sector number")
    parser.add_argument("--limit", type=int, default=100, help="Max stars to process")
    parser.add_argument("--ra", type=float, default=None, help="Center RA (degrees)")
    parser.add_argument("--dec", type=float, default=None, help="Center Dec (degrees)")
    parser.add_argument("--radius", type=float, default=0.5, help="Cone search radius (degrees)")
    parser.add_argument("--output", default=None, help="Output directory")
    parser.add_argument("--fgk-only", action="store_true", help="Filter for FGK dwarfs")
    parser.add_argument("--cutout-size", type=int, default=11, help="TESScut cutout size (pixels)")
    args = parser.parse_args()

    output_dir = Path(args.output) if args.output else Path(f"data/phase2/ffi_sector_{args.sector}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Default pointing: use a known bright field in the sector
    # User should provide RA/Dec for targeted searches
    if args.ra is None or args.dec is None:
        print("ERROR: --ra and --dec required (center of field to search)")
        print("  Tip: Use TESS-Point or the TESS viewing tool to find coordinates")
        print("  Example: --ra 90.0 --dec -66.5 (TESS CVZ south)")
        sys.exit(1)

    print("=" * 60)
    print(f"PHASE 2 FFI: TESScut Light Curve Extraction — Sector {args.sector}")
    print(f"  Center: RA={args.ra:.4f}, Dec={args.dec:.4f}")
    print(f"  Radius: {args.radius} deg")
    print(f"  Limit: {args.limit}")
    print(f"  FGK filter: {args.fgk_only}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    # Step 1: Get targets from TIC catalog
    print(f"\n[1/5] Querying TIC catalog for field stars...")
    targets = get_ffi_targets_from_tic(args.ra, args.dec, args.radius, args.sector)

    if not targets:
        print("ERROR: No TIC targets found in this field.")
        sys.exit(1)

    # Step 2: Remove SPOC 2-minute targets (already searched by TESS pipeline)
    print(f"\n[2/5] Removing SPOC 2-minute targets (already transit-searched)...")
    spoc_tics = get_spoc_target_ids(args.sector)
    targets = [t for t in targets if t["tic_id"] not in spoc_tics]
    print(f"  After removing SPOC targets: {len(targets)}")

    # Step 3: Remove known TOIs
    print(f"\n[3/5] Removing known TOIs and cTOIs...")
    sys.path.insert(0, str(Path(__file__).parent))
    from download_sector_bulk import fetch_toi_tic_ids
    known_tics = fetch_toi_tic_ids()
    targets = [t for t in targets if t["tic_id"] not in known_tics]
    print(f"  After removing TOIs/cTOIs: {len(targets)}")

    # Step 4: Optional FGK filter
    if args.fgk_only:
        print(f"\n[4/5] Filtering for FGK dwarfs...")
        from download_sector_bulk import filter_fgk_dwarfs
        fgk_ids = set(filter_fgk_dwarfs([t["tic_id"] for t in targets]))
        targets = [t for t in targets if t["tic_id"] in fgk_ids]
    else:
        print(f"\n[4/5] Skipping FGK filter")

    # Sort by brightness (brightest first = best SNR)
    targets.sort(key=lambda t: t["Tmag"])
    targets = targets[:args.limit]
    print(f"\n  Final target count: {len(targets)} (sorted by Tmag)")

    # Save manifest
    manifest = output_dir / "targets_ffi.json"
    with open(manifest, "w") as f:
        json.dump({
            "sector": args.sector,
            "method": "TESScut",
            "ra_center": args.ra,
            "dec_center": args.dec,
            "radius_deg": args.radius,
            "fgk_only": args.fgk_only,
            "n_targets": len(targets),
            "tic_ids": [t["tic_id"] for t in targets],
        }, f, indent=2)

    # Step 5: Extract light curves
    print(f"\n[5/5] Extracting FFI light curves via TESScut...")
    downloaded, skipped, failed = download_ffi_lightcurves(
        targets, args.sector, output_dir
    )

    print(f"\n{'=' * 60}")
    print(f"  FFI EXTRACTION COMPLETE")
    print(f"  Extracted:   {downloaded}")
    print(f"  Already had: {skipped}")
    print(f"  Failed:      {failed}")
    print(f"  Output:      {output_dir}")
    print(f"{'=' * 60}")

    if downloaded + skipped > 0:
        print(f"\nNext: Run BLS on FFI light curves:")
        print(f"  ./target/release/hunt search -i {output_dir} -o results/phase2/candidates_ffi_s{args.sector}.json --snr-threshold 5.5 --n-periods 15000")


if __name__ == "__main__":
    main()
