#!/usr/bin/env python3.11
"""Phase 2: Flag potential new planet discoveries from BLS candidates.

Cross-references BLS candidates against:
  1. ExoFOP TOI catalog (known TESS Objects of Interest)
  2. ExoFOP community TOI catalog (cTOIs)
  3. NASA Confirmed Exoplanet catalog
  4. Validation results (planet_score threshold)

Outputs a ranked discovery report with candidates categorized as:
  - NEW: Not in any catalog, passes validation — submit to ExoFOP as cTOI
  - KNOWN_TOI: Already a TOI — skip
  - KNOWN_PLANET: Already confirmed — skip
  - LOW_SCORE: Below validation threshold — likely false positive

Usage:
  python3.11 python/flag_discoveries.py \
    --input results/phase2/candidates_s56.json \
    --validation results/phase2/validation_results.json \
    --min-score 60
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")


def extract_tic_id(filename: str) -> int | None:
    """Extract TIC ID from a filename like 'TIC_261136679_s0056.csv'."""
    match = re.search(r"TIC[_\s]?(\d+)", filename, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def fetch_known_catalogs() -> tuple[set[int], set[int], set[int]]:
    """Fetch TOI, cTOI, and confirmed planet TIC IDs.

    Returns (toi_tics, ctoi_tics, confirmed_tics).
    """
    toi_tics = set()
    ctoi_tics = set()
    confirmed_tics = set()

    # TOIs
    try:
        toi_url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
        tois = pd.read_csv(toi_url, comment="#")
        toi_tics = set(tois["TIC ID"].dropna().astype(int).unique())
        print(f"  TOI catalog: {len(toi_tics)} unique TIC IDs")
    except Exception as e:
        print(f"  WARNING: Could not fetch TOI catalog: {e}")

    # Community TOIs
    try:
        ctoi_url = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"
        ctois = pd.read_csv(ctoi_url, comment="#")
        ctoi_tics = set(ctois["TIC ID"].dropna().astype(int).unique())
        print(f"  cTOI catalog: {len(ctoi_tics)} unique TIC IDs")
    except Exception as e:
        print(f"  WARNING: Could not fetch cTOI catalog: {e}")

    # Confirmed planets (from local file if available)
    local_catalog = Path("data/lightcurves/confirmed_exoplanets.csv")
    if local_catalog.exists():
        try:
            confirmed = pd.read_csv(local_catalog)
            if "tic_id" in confirmed.columns:
                confirmed_tics = set(confirmed["tic_id"].dropna().astype(int).unique())
            print(f"  Confirmed catalog: {len(confirmed_tics)} TIC IDs")
        except Exception:
            pass

    return toi_tics, ctoi_tics, confirmed_tics


def classify_candidate(
    tic_id: int | None,
    planet_score: int | None,
    min_score: int,
    toi_tics: set[int],
    ctoi_tics: set[int],
    confirmed_tics: set[int],
) -> str:
    """Classify a candidate as NEW, KNOWN_TOI, KNOWN_PLANET, or LOW_SCORE."""
    if tic_id is not None:
        if tic_id in confirmed_tics:
            return "KNOWN_PLANET"
        if tic_id in toi_tics or tic_id in ctoi_tics:
            return "KNOWN_TOI"

    if planet_score is not None and planet_score < min_score:
        return "LOW_SCORE"

    return "NEW"


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2: Flag new planet discoveries from BLS candidates"
    )
    parser.add_argument(
        "--input", required=True,
        help="BLS candidates JSON from hunt search"
    )
    parser.add_argument(
        "--validation", default=None,
        help="Validation results JSON from hunt validate"
    )
    parser.add_argument(
        "--min-score", type=int, default=60,
        help="Minimum planet_score to consider (default: 60)"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output report path (default: same dir as input, discovery_report.json)"
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: {input_path} not found")
        sys.exit(1)

    print("=" * 60)
    print("PHASE 2: Discovery Flagging")
    print("=" * 60)

    # Load BLS candidates
    with open(input_path) as f:
        hunt_report = json.load(f)

    candidates = hunt_report.get("candidates", [])
    print(f"\nLoaded {len(candidates)} BLS candidates from {input_path}")

    # Load validation results if available
    scores = {}
    if args.validation:
        val_path = Path(args.validation)
        if val_path.exists():
            with open(val_path) as f:
                val_results = json.load(f)
            for v in val_results:
                scores[v["filename"]] = v.get("planet_score", 0)
            print(f"Loaded {len(scores)} validation scores from {val_path}")

    # Fetch known catalogs
    print("\nFetching known catalogs...")
    toi_tics, ctoi_tics, confirmed_tics = fetch_known_catalogs()
    all_known = toi_tics | ctoi_tics | confirmed_tics
    print(f"  Total known TIC IDs: {len(all_known)}")

    # Classify each candidate
    print(f"\nClassifying candidates (min_score={args.min_score})...")
    results = []
    for c in candidates:
        tic_id = extract_tic_id(c["filename"])
        planet_score = scores.get(c["filename"])
        status = classify_candidate(
            tic_id, planet_score, args.min_score,
            toi_tics, ctoi_tics, confirmed_tics,
        )
        results.append({
            "filename": c["filename"],
            "tic_id": tic_id,
            "period_days": c["period_days"],
            "depth_ppm": c["depth_ppm"],
            "snr": c["snr"],
            "radius_ratio": c["radius_ratio"],
            "n_transits": c["n_transits"],
            "planet_score": planet_score,
            "status": status,
        })

    # Summary
    new = [r for r in results if r["status"] == "NEW"]
    known_toi = [r for r in results if r["status"] == "KNOWN_TOI"]
    known_planet = [r for r in results if r["status"] == "KNOWN_PLANET"]
    low_score = [r for r in results if r["status"] == "LOW_SCORE"]

    print(f"\n{'=' * 60}")
    print(f"  DISCOVERY REPORT")
    print(f"{'=' * 60}")
    print(f"  Total candidates:   {len(results)}")
    print(f"  NEW (potential):    {len(new)}")
    print(f"  KNOWN_TOI:          {len(known_toi)}")
    print(f"  KNOWN_PLANET:       {len(known_planet)}")
    print(f"  LOW_SCORE:          {len(low_score)}")
    print(f"{'=' * 60}")

    if new:
        # Sort by planet_score (descending), then SNR
        new.sort(key=lambda r: (-(r["planet_score"] or 0), -r["snr"]))
        print(f"\n  TOP NEW CANDIDATES:")
        print(f"  {'TIC ID':>12}  {'Period(d)':>10}  {'SNR':>8}  {'Depth(ppm)':>11}  {'Rp/Rs':>8}  {'Score':>6}")
        print(f"  {'─' * 62}")
        for r in new[:20]:
            tic = r["tic_id"] or "?"
            score = r["planet_score"] if r["planet_score"] is not None else "?"
            print(f"  {tic:>12}  {r['period_days']:>10.4f}  {r['snr']:>8.1f}  {r['depth_ppm']:>11.0f}  {r['radius_ratio']:>8.4f}  {score:>6}")

        # Flag planet-sized candidates specifically
        planet_sized = [r for r in new if r["radius_ratio"] < 0.3]
        if planet_sized:
            print(f"\n  PLANET-SIZED NEW CANDIDATES (Rp/Rs < 0.3):")
            for r in planet_sized:
                tic = r["tic_id"] or "?"
                rp_earth = r["radius_ratio"] * 1.0  # approximate, need stellar radius
                print(f"    TIC {tic}: P={r['period_days']:.4f}d, SNR={r['snr']:.1f}, Rp/Rs={r['radius_ratio']:.4f}")
            print(f"\n  These are your best candidates for ExoFOP cTOI submission!")
        else:
            print(f"\n  No planet-sized (Rp/Rs < 0.3) new candidates found.")
            print(f"  All detections have Rp/Rs > 0.3, likely eclipsing binaries.")
    else:
        print("\n  No new candidates found in this batch.")

    # Save report
    output_path = Path(args.output) if args.output else input_path.parent / "discovery_report.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump({
            "summary": {
                "total": len(results),
                "new": len(new),
                "known_toi": len(known_toi),
                "known_planet": len(known_planet),
                "low_score": len(low_score),
                "min_score_threshold": args.min_score,
            },
            "candidates": results,
        }, f, indent=2)
    print(f"\n  Report saved: {output_path}")

    if new:
        print(f"\n  Next steps for NEW candidates:")
        print(f"  1. Generate phase-fold plots for visual inspection")
        print(f"  2. Run deep validation (TLS, TRICERATOPS) on top candidates")
        print(f"  3. Check MAST/Simbad for known variables at these coordinates")
        print(f"  4. Submit as cTOI to ExoFOP: https://exofop.ipac.caltech.edu/tess/")


if __name__ == "__main__":
    main()
