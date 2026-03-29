#!/usr/bin/env python3
"""
🔭 download_lightcurves.py — Fetch TESS/Kepler light curves for transit hunting.

This script downloads light curves from MAST (Mikulski Archive for Space Telescopes)
and converts them to simple CSV files that the Rust BLS engine can process.

Usage:
    python download_lightcurves.py --sector 56 --limit 500
    python download_lightcurves.py --target "TIC 261136679" --mission tess
    python download_lightcurves.py --candidates-only  # download only TOIs (unconfirmed)

Requirements:
    pip install lightkurve astroquery pandas numpy tqdm
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


def setup_dirs(output_dir: str) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out


def download_tess_sector(sector: int, output_dir: Path, limit: int = 200, cadence: str = "short") -> int:
    """Download TESS light curves for a given sector."""
    import lightkurve as lk

    print(f"\n🛰️  Searching TESS Sector {sector} (limit={limit})...")
    search = lk.search_lightcurve(
        f"sector {sector}",
        mission="TESS",
        author="SPOC",  # Science Processing Operations Center pipeline
        cadence=cadence,
    )

    if len(search) == 0:
        print(f"   ⚠️  No results for sector {sector}. Try a different sector.")
        return 0

    print(f"   Found {len(search)} light curves. Downloading up to {limit}...")

    downloaded = 0
    for i, result in enumerate(tqdm(search[:limit], desc="   Downloading")):
        try:
            lc = result.download(quality_bitmask="hardest")
            if lc is None:
                continue

            # Normalize and clean
            lc = lc.remove_nans().remove_outliers(sigma=5)
            lc = lc.normalize()

            if len(lc.time.value) < 100:
                continue

            # Save as CSV
            target_name = str(result.target_name).replace(" ", "_").replace("/", "_")
            filename = f"{target_name}_s{sector:04d}.csv"
            filepath = output_dir / filename

            df = pd.DataFrame({
                "time": lc.time.value,
                "flux": lc.flux.value,
                "flux_err": lc.flux_err.value if lc.flux_err is not None else 0.001,
            })
            df.to_csv(filepath, index=False)
            downloaded += 1

        except (OSError, ConnectionError, ValueError, RuntimeError) as e:
            print(f"   Warning: skipping {result.target_name}: {e}")
            continue

    print(f"   ✅ Downloaded {downloaded} light curves to {output_dir}/")
    return downloaded


def download_toi_candidates(output_dir: Path, limit: int = 500) -> int:
    """Download light curves for TESS Objects of Interest (unconfirmed candidates).

    These are the most interesting targets — planets waiting to be confirmed!
    """
    import lightkurve as lk

    print("\n🎯 Fetching TESS Objects of Interest (TOIs)...")

    # Get TOI list from ExoFOP
    try:
        toi_url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        tois = pd.read_csv(toi_url, comment="#")
        print(f"   Found {len(tois)} TOIs in catalog")

        # Filter for interesting candidates
        # PC = Planet Candidate, KP = Known Planet, CP = Confirmed Planet
        candidates = tois[tois["TFOPWG Disposition"].isin(["PC", ""])]
        candidates = candidates.head(limit)

        print(f"   Downloading light curves for {len(candidates)} unconfirmed candidates...")

    except Exception as e:
        print(f"   ⚠️  Could not fetch TOI list: {e}")
        print("   Falling back to sector-based download...")
        return download_tess_sector(1, output_dir, limit)

    downloaded = 0
    for _, row in tqdm(candidates.iterrows(), total=len(candidates), desc="   Downloading"):
        try:
            tic_id = f"TIC {int(row['TIC ID'])}"
            search = lk.search_lightcurve(tic_id, mission="TESS", author="SPOC")

            if len(search) == 0:
                continue

            lc = search[0].download(quality_bitmask="hardest")
            if lc is None:
                continue

            lc = lc.remove_nans().remove_outliers(sigma=5).normalize()

            if len(lc.time.value) < 100:
                continue

            target_name = tic_id.replace(" ", "_")
            toi_num = str(row.get("TOI", "unknown"))
            filename = f"TOI_{toi_num}_{target_name}.csv"
            filepath = output_dir / filename

            df = pd.DataFrame({
                "time": lc.time.value,
                "flux": lc.flux.value,
                "flux_err": lc.flux_err.value if lc.flux_err is not None else 0.001,
            })
            df.to_csv(filepath, index=False)
            downloaded += 1

        except (OSError, ConnectionError, ValueError, RuntimeError) as e:
            print(f"   Warning: skipping TIC {int(row.get('TIC ID', 0))}: {e}")
            continue

    print(f"   ✅ Downloaded {downloaded} TOI light curves to {output_dir}/")
    return downloaded


def download_kepler_kois(output_dir: Path, limit: int = 200) -> int:
    """Download Kepler Objects of Interest — the original planet hunting dataset."""
    import lightkurve as lk

    print("\n🔭 Fetching Kepler Objects of Interest (KOIs)...")

    # Get the KOI table from NASA Exoplanet Archive
    try:
        from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

        kois = NasaExoplanetArchive.query_criteria(
            table="koi",
            select="kepoi_name,kepid,koi_disposition,koi_period,koi_depth",
            where="koi_disposition='CANDIDATE'",
            order="koi_depth desc",
        )
        koi_ids = [f"KIC {kid}" for kid in kois["kepid"][:limit]]
        print(f"   Found {len(koi_ids)} KOI candidates")
    except (OSError, ConnectionError) as e:
        print("   Using fallback KOI list...")
        koi_ids = [f"KIC {kid}" for kid in [8191672, 3558849, 5728139, 10797460, 7040629]]

    downloaded = 0
    for kic in tqdm(koi_ids[:limit], desc="   Downloading"):
        try:
            search = lk.search_lightcurve(kic, mission="Kepler", author="Kepler")
            if len(search) == 0:
                continue

            lc = search[0].download(quality_bitmask="hardest")
            if lc is None:
                continue

            lc = lc.remove_nans().remove_outliers(sigma=5).normalize()

            if len(lc.time.value) < 100:
                continue

            target_name = kic.replace(" ", "_")
            filename = f"{target_name}_kepler.csv"
            filepath = output_dir / filename

            df = pd.DataFrame({
                "time": lc.time.value,
                "flux": lc.flux.value,
                "flux_err": lc.flux_err.value if lc.flux_err is not None else 0.001,
            })
            df.to_csv(filepath, index=False)
            downloaded += 1

        except (OSError, ConnectionError, ValueError, RuntimeError) as e:
            print(f"   Warning: skipping {kic}: {e}")
            continue

    print(f"   ✅ Downloaded {downloaded} KOI light curves to {output_dir}/")
    return downloaded


def download_exoplanet_catalog(output_dir: Path) -> None:
    """Download the full confirmed exoplanet catalog from NASA for cross-referencing."""
    print("\n📊 Downloading NASA confirmed exoplanet catalog...")

    try:
        from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

        planets = NasaExoplanetArchive.query_criteria(
            table="pscomppars",
            select="pl_name,hostname,pl_orbper,pl_rade,pl_bmasse,pl_eqt,st_teff,st_rad,sy_dist,disc_facility",
        )
        df = planets.to_pandas()
        filepath = output_dir / "confirmed_exoplanets.csv"
        df.to_csv(filepath, index=False)
        print(f"   ✅ Saved {len(df)} confirmed exoplanets to {filepath}")

    except Exception as e:
        print(f"   ⚠️  Could not download catalog: {e}")
        print("   You can manually download from: https://exoplanetarchive.ipac.caltech.edu/")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="🔭 Download light curves for exoplanet hunting"
    )
    parser.add_argument(
        "--mission", choices=["tess", "kepler", "both"], default="tess",
        help="Which mission's data to download"
    )
    parser.add_argument("--sector", type=int, default=56, help="TESS sector number (1-72+)")
    parser.add_argument("--limit", type=int, default=200, help="Max light curves to download")
    parser.add_argument("--output", default="data/lightcurves", help="Output directory")
    parser.add_argument(
        "--candidates-only", action="store_true",
        help="Only download unconfirmed TOIs/KOIs (best for discovery!)"
    )
    parser.add_argument(
        "--catalog", action="store_true",
        help="Also download the confirmed exoplanet catalog"
    )

    args = parser.parse_args()
    output_dir = setup_dirs(args.output)

    print("🔭 Exoplanet Hunter — Data Downloader")
    print("━" * 50)

    total = 0

    if args.candidates_only:
        if args.mission in ("tess", "both"):
            total += download_toi_candidates(output_dir, args.limit)
        if args.mission in ("kepler", "both"):
            total += download_kepler_kois(output_dir, args.limit)
    else:
        if args.mission in ("tess", "both"):
            total += download_tess_sector(args.sector, output_dir, args.limit)
        if args.mission in ("kepler", "both"):
            total += download_kepler_kois(output_dir, args.limit)

    if args.catalog:
        download_exoplanet_catalog(output_dir)

    print(f"\n{'━' * 50}")
    print(f"🏁 Total: {total} light curves downloaded")
    print(f"📂 Location: {output_dir.resolve()}")
    print(f"\n🚀 Next step: Run the Rust BLS detector:")
    print(f"   cargo run --release -- -i {output_dir} -o candidates.json")


if __name__ == "__main__":
    main()
