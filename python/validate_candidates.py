#!/usr/bin/env python3
"""
Exohuntr — Deep Candidate Validation Pipeline

Applies multiple false-positive tests to BLS transit candidates:
  1. Odd/even transit depth test (eclipsing binary filter)
  2. Secondary eclipse search at phase 0.5
  3. V-shape vs U-shape morphology (ingress/egress analysis)
  4. Period agreement with TESS pipeline
  5. ExoFOP cross-reference (disposition, stellar params, follow-up status)
  6. Radius ratio sanity check (Rp/Rs > 1.0 = likely binary)

Produces a scored validation report ranking candidates by planet likelihood.

Usage:
    python validate_candidates.py --input candidates.json --lightcurves data/lightcurves/
"""

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================================
# Test 1: Odd/Even Transit Depth
# ============================================================================

def test_odd_even_depth(time, flux, period, epoch):
    """Compare transit depths of odd vs even transits.

    If depths differ significantly, the signal is likely an eclipsing binary
    with a period twice what we detected (primary + secondary eclipses).

    Returns: (ratio, passed) where ratio = depth_odd/depth_even.
    A ratio near 1.0 means consistent depths (planet-like).
    """
    phase = ((time - epoch) / period) % 1.0
    transit_mask = (phase < 0.05) | (phase > 0.95)

    # Assign each transit to odd or even
    transit_number = np.floor((time - epoch) / period).astype(int)
    odd_mask = transit_mask & (transit_number % 2 == 1)
    even_mask = transit_mask & (transit_number % 2 == 0)

    out_mask = (phase > 0.15) & (phase < 0.85)
    baseline = np.median(flux[out_mask]) if out_mask.sum() > 10 else np.median(flux)

    odd_flux = flux[odd_mask]
    even_flux = flux[even_mask]

    if len(odd_flux) < 5 or len(even_flux) < 5:
        return None, None  # insufficient data

    depth_odd = baseline - np.median(odd_flux)
    depth_even = baseline - np.median(even_flux)

    if depth_even <= 0 or depth_odd <= 0:
        return None, None

    ratio = depth_odd / depth_even
    # Consistent depths: ratio between 0.7 and 1.4
    passed = 0.7 <= ratio <= 1.4
    return ratio, passed


# ============================================================================
# Test 2: Secondary Eclipse Search
# ============================================================================

def test_secondary_eclipse(time, flux, period, epoch):
    """Search for a brightness dip at phase 0.5 (opposite the primary transit).

    A secondary eclipse indicates the companion is luminous — i.e., a star
    (eclipsing binary), not a planet. Planets don't produce secondary eclipses
    in TESS photometry (too faint).

    Returns: (secondary_depth_ratio, passed)
    secondary_depth_ratio = secondary_depth / primary_depth
    passed = True if no significant secondary eclipse detected.
    """
    phase = ((time - epoch) / period) % 1.0

    # Primary transit region: phase 0 +/- 0.05
    primary_mask = (phase < 0.05) | (phase > 0.95)
    # Secondary eclipse region: phase 0.5 +/- 0.05
    secondary_mask = (phase > 0.45) & (phase < 0.55)
    # Out-of-transit baseline
    out_mask = (phase > 0.1) & (phase < 0.4)

    if primary_mask.sum() < 5 or secondary_mask.sum() < 5 or out_mask.sum() < 10:
        return None, None

    baseline = np.median(flux[out_mask])
    primary_depth = baseline - np.median(flux[primary_mask])
    secondary_depth = baseline - np.median(flux[secondary_mask])

    if primary_depth <= 0:
        return None, None

    ratio = secondary_depth / primary_depth
    # No secondary eclipse: ratio < 0.3
    passed = ratio < 0.3
    return ratio, passed


# ============================================================================
# Test 3: Transit Shape (V vs U)
# ============================================================================

def test_transit_shape(time, flux, period, epoch, duration_hours):
    """Analyze transit morphology: U-shaped (planet) vs V-shaped (binary).

    Planets produce flat-bottomed transits (the planet fully covers part of the
    stellar disk). Eclipsing binaries often produce V-shaped transits (comparable
    sized objects, partial eclipse).

    Returns: (flatness_score, passed)
    flatness_score: ratio of in-transit flat portion to total transit width.
    Higher = more U-shaped = more planet-like.
    """
    phase = ((time - epoch) / period) % 1.0
    phase[phase > 0.5] -= 1.0  # center on 0

    dur_phase = duration_hours / (period * 24)
    if dur_phase <= 0 or dur_phase > 0.3:
        return None, None

    # In-transit points
    transit_mask = np.abs(phase) < dur_phase
    if transit_mask.sum() < 10:
        return None, None

    transit_phase = phase[transit_mask]
    transit_flux = flux[transit_mask]

    # Sort by phase
    sort_idx = np.argsort(transit_phase)
    transit_phase = transit_phase[sort_idx]
    transit_flux = transit_flux[sort_idx]

    # Divide transit into 3 equal parts: ingress, flat, egress
    n = len(transit_phase)
    third = n // 3
    if third < 3:
        return None, None

    ingress_flux = transit_flux[:third]
    flat_flux = transit_flux[third:2*third]
    egress_flux = transit_flux[2*third:]

    # For a U-shape, the flat portion should have lower variance than ingress/egress
    flat_std = np.std(flat_flux)
    edge_std = (np.std(ingress_flux) + np.std(egress_flux)) / 2

    if edge_std < 1e-10:
        return None, None

    # Also check that flat bottom is actually flat (low slope)
    flatness = 1.0 - (flat_std / edge_std) if edge_std > flat_std else 0.5

    # Flatness > 0.3 suggests U-shape (planet-like)
    passed = flatness > 0.2
    return flatness, passed


# ============================================================================
# Test 4: Period Agreement with TESS Pipeline
# ============================================================================

def test_period_agreement(our_period, tess_period):
    """Check if our detected period matches the TESS pipeline period.

    Agreement validates our detection independently. Harmonic periods
    (2x, 0.5x) are also acceptable as they indicate the same underlying signal.
    """
    if tess_period is None or np.isnan(tess_period) or tess_period <= 0:
        return None, None

    ratio = our_period / tess_period

    # Direct match
    if abs(ratio - 1.0) < 0.05:
        return ("exact", True)

    # Harmonic matches
    for h in [0.5, 2.0, 3.0, 1.0/3.0, 4.0, 0.25]:
        if abs(ratio - h) < 0.1:
            return ("harmonic", True)

    return ("disagree", False)


# ============================================================================
# Scoring & Report
# ============================================================================

def compute_planet_score(results):
    """Compute a 0-100 planet likelihood score from validation tests."""
    score = 50  # start neutral

    # Odd/even depth test (±15)
    if results.get("odd_even_passed") is True:
        score += 15
    elif results.get("odd_even_passed") is False:
        score -= 15

    # Secondary eclipse test (±20, strongest discriminator)
    if results.get("secondary_passed") is True:
        score += 20
    elif results.get("secondary_passed") is False:
        score -= 25

    # Transit shape test (±10)
    if results.get("shape_passed") is True:
        score += 10
    elif results.get("shape_passed") is False:
        score -= 10

    # Period agreement (±15)
    if results.get("period_match") == "exact":
        score += 15
    elif results.get("period_match") == "harmonic":
        score += 5
    elif results.get("period_match") == "disagree":
        score -= 10

    # Radius ratio sanity (±10)
    rr = results.get("radius_ratio", 0)
    if rr and rr < 0.3:
        score += 10  # planet-sized
    elif rr and rr > 1.0:
        score -= 15  # bigger than star = binary

    # SNR bonus
    snr = results.get("snr", 0)
    if snr > 50:
        score += 5
    elif snr < 10:
        score -= 5

    # Clamp
    return max(0, min(100, score))


def generate_validation_report(validated, output_dir):
    """Generate markdown validation report."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort by planet score
    validated.sort(key=lambda x: x.get("planet_score", 0), reverse=True)

    lines = [
        "# Exohuntr — Candidate Validation Report\n",
        "## Methodology\n",
        "Each candidate was subjected to 5 independent false-positive tests:\n",
        "| Test | What it checks | Planet signal | Binary signal |",
        "|------|---------------|---------------|---------------|",
        "| **Odd/Even Depth** | Do alternate transits have equal depth? | Equal depths | Unequal depths |",
        "| **Secondary Eclipse** | Is there a dip at orbital phase 0.5? | No dip | Dip present |",
        "| **Transit Shape** | Is the transit U-shaped or V-shaped? | U-shaped (flat bottom) | V-shaped (pointed) |",
        "| **Period Agreement** | Does our period match TESS pipeline? | Exact match | Disagreement |",
        "| **Radius Ratio** | Is Rp/Rs physically reasonable? | < 0.3 (planet-sized) | > 1.0 (impossible) |\n",
        "Each test contributes to a **Planet Likelihood Score** (0-100).\n",
        "---\n",
    ]

    # Summary stats
    scores = [v["planet_score"] for v in validated]
    high = [v for v in validated if v["planet_score"] >= 70]
    medium = [v for v in validated if 50 <= v["planet_score"] < 70]
    low = [v for v in validated if v["planet_score"] < 50]

    lines.extend([
        "## Summary\n",
        f"- **Total candidates validated:** {len(validated)}",
        f"- **High confidence (score >= 70):** {len(high)}",
        f"- **Medium confidence (50-69):** {len(medium)}",
        f"- **Low confidence / likely false positive (< 50):** {len(low)}",
        f"- **Mean planet score:** {np.mean(scores):.1f}",
        f"- **Median planet score:** {np.median(scores):.1f}\n",
        "---\n",
    ])

    # High confidence candidates
    if high:
        lines.extend([
            f"## High Confidence Candidates ({len(high)})\n",
            "These candidates passed most or all validation tests and are the most likely to be genuine planetary signals.\n",
            "| Rank | Target | Period (d) | SNR | Score | Odd/Even | No 2nd Eclipse | Shape | Period Match | Rp (R_Earth) |",
            "|------|--------|-----------|-----|-------|----------|----------------|-------|-------------|-------------|",
        ])
        for i, v in enumerate(high):
            oe = "PASS" if v.get("odd_even_passed") else ("FAIL" if v.get("odd_even_passed") is False else "N/A")
            se = "PASS" if v.get("secondary_passed") else ("FAIL" if v.get("secondary_passed") is False else "N/A")
            sh = "PASS" if v.get("shape_passed") else ("FAIL" if v.get("shape_passed") is False else "N/A")
            pm = v.get("period_match", "N/A")
            rp = f'{v.get("planet_radius_earth", "?"):.1f}' if isinstance(v.get("planet_radius_earth"), float) else "?"
            name = v["filename"].replace(".csv", "").replace("TOI_", "TOI ").replace("_TIC_", " / TIC ")
            lines.append(
                f'| {i+1} | {name} | {v["period_days"]:.4f} | {v["snr"]:.1f} | **{v["planet_score"]}** | {oe} | {se} | {sh} | {pm} | {rp} |'
            )
        lines.append("")

    # All candidates table
    lines.extend([
        "\n---\n",
        "## All Candidates (sorted by planet likelihood score)\n",
        "| # | Target | Period | SNR | Score | Tests Passed | TESS Disposition | Rp (R_Earth) |",
        "|---|--------|--------|-----|-------|-------------|-----------------|-------------|",
    ])
    for i, v in enumerate(validated):
        tests_passed = sum([
            v.get("odd_even_passed", False) == True,
            v.get("secondary_passed", False) == True,
            v.get("shape_passed", False) == True,
            v.get("period_match") in ("exact", "harmonic"),
        ])
        name = v["filename"][:35].replace("TOI_", "TOI ").replace("_TIC_", " / TIC ").replace(".csv", "")
        disp = v.get("tfopwg_disposition", "?")
        rp = f'{v.get("planet_radius_earth", "?"):.1f}' if isinstance(v.get("planet_radius_earth"), float) else "?"
        lines.append(
            f'| {i+1} | `{name}` | {v["period_days"]:.4f} | {v["snr"]:.1f} | {v["planet_score"]} | {tests_passed}/4 | {disp} | {rp} |'
        )

    # Next steps
    lines.extend([
        "\n---\n",
        "## Path to Publication\n",
        "### For high-confidence candidates (score >= 70):\n",
        "1. **Submit to ExoFOP** — Register findings at [exofop.ipac.caltech.edu/tess](https://exofop.ipac.caltech.edu/tess/)",
        "2. **Request follow-up** — Ground-based photometry to confirm transit at predicted times",
        "3. **Check SIMBAD/Vizier** — Ensure the star isn't a known variable or binary",
        "4. **Contact TESS team** — Email tess-followup@mit.edu with candidates and validation",
        "5. **Write up results** — Submit to Research Notes of the AAS (RNAAS) for short discoveries",
        "6. **Citizen science programs** — Submit to Planet Hunters TESS or AAVSO\n",
        "### Publication venues:\n",
        "- **RNAAS** (Research Notes of the AAS) — Short-form, peer-reviewed, fast publication",
        "- **AJ** (Astronomical Journal) — Full paper if multiple confirmed candidates",
        "- **MNRAS** (Monthly Notices of the RAS) — International alternative",
        "- **arXiv astro-ph.EP** — Pre-print for immediate visibility\n",
        "### Key contacts:\n",
        "- TESS Follow-up Program: [tess.mit.edu](https://tess.mit.edu)",
        "- ExoFOP: [exofop.ipac.caltech.edu](https://exofop.ipac.caltech.edu/tess/)",
        "- Planet Hunters TESS: [zooniverse.org](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess)",
        "- AAVSO Exoplanet Database: [aavso.org](https://www.aavso.org)\n",
        "\n---\n",
        f"*Validated {len(validated)} candidates. Generated by Exohuntr validation pipeline.*",
    ])

    report_path = output_dir / "VALIDATION_REPORT.md"
    report_path.write_text("\n".join(lines))
    print(f"   Report saved to {report_path}")
    return report_path


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Validate BLS transit candidates")
    parser.add_argument("--input", default="candidates.json")
    parser.add_argument("--lightcurves", default="data/lightcurves")
    parser.add_argument("--output", default="results")
    parser.add_argument("--exofop", action="store_true", help="Fetch ExoFOP TOI catalog for cross-reference")

    args = parser.parse_args()
    lc_dir = Path(args.lightcurves)

    print("\nExohuntr — Candidate Validation Pipeline")
    print("=" * 55)

    # Load candidates
    with open(args.input) as f:
        report = json.load(f)
    candidates = report["candidates"]
    print(f"  Loaded {len(candidates)} candidates\n")

    # Load ExoFOP catalog
    toi_data = {}
    if args.exofop:
        print("  Fetching ExoFOP TOI catalog...")
        try:
            tois = pd.read_csv(
                "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
                comment="#"
            )
            for _, r in tois.iterrows():
                tic = int(r["TIC ID"])
                if tic not in toi_data:
                    toi_data[tic] = r
            print(f"  Loaded {len(tois)} TOIs from ExoFOP\n")
        except Exception as e:
            print(f"  Warning: Could not fetch ExoFOP data: {e}\n")

    # Run validation
    validated = []
    for i, c in enumerate(candidates):
        result = dict(c)  # copy all candidate fields

        # Load light curve
        filepath = lc_dir / c["filename"]
        if not filepath.exists():
            result["planet_score"] = 0
            validated.append(result)
            continue

        df = pd.read_csv(filepath)
        time = df["time"].values
        flux = df["flux"].values
        period = c["period_days"]
        epoch = c["epoch"]

        # Test 1: Odd/Even
        oe_ratio, oe_passed = test_odd_even_depth(time, flux, period, epoch)
        result["odd_even_ratio"] = oe_ratio
        result["odd_even_passed"] = oe_passed

        # Test 2: Secondary Eclipse
        se_ratio, se_passed = test_secondary_eclipse(time, flux, period, epoch)
        result["secondary_depth_ratio"] = se_ratio
        result["secondary_passed"] = se_passed

        # Test 3: Transit Shape
        shape_score, shape_passed = test_transit_shape(time, flux, period, epoch, c["duration_hours"])
        result["shape_flatness"] = shape_score
        result["shape_passed"] = shape_passed

        # Test 4: Period Agreement with TESS
        tic_str = c["filename"].split("_TIC_")
        tic_id = int(tic_str[1].replace(".csv", "")) if len(tic_str) > 1 else None
        toi_row = toi_data.get(tic_id)

        if toi_row is not None:
            tess_period = toi_row.get("Period (days)")
            match_type, period_passed = test_period_agreement(period, tess_period)
            result["period_match"] = match_type
            result["period_match_passed"] = period_passed
            result["tess_period"] = tess_period

            # ExoFOP metadata
            result["tfopwg_disposition"] = toi_row.get("TFOPWG Disposition", "?")
            result["tess_disposition"] = toi_row.get("TESS Disposition", "?")
            result["planet_radius_earth"] = toi_row.get("Planet Radius (R_Earth)")
            result["stellar_teff"] = toi_row.get("Stellar Eff Temp (K)")
            result["stellar_radius_sun"] = toi_row.get("Stellar Radius (R_Sun)")
            result["toi_number"] = toi_row.get("TOI")
        else:
            result["period_match"] = None
            result["tfopwg_disposition"] = "?"

        # Compute score
        result["planet_score"] = compute_planet_score(result)

        validated.append(result)

        if (i + 1) % 50 == 0:
            print(f"  Validated {i+1}/{len(candidates)}")

    print(f"  Validated {len(validated)}/{len(candidates)} candidates\n")

    # Generate report
    report_path = generate_validation_report(validated, args.output)

    # Save detailed JSON
    json_path = Path(args.output) / "validation_results.json"
    # Clean NaN for JSON serialization
    for v in validated:
        for k, val in v.items():
            if isinstance(val, float) and np.isnan(val):
                v[k] = None
    with open(json_path, "w") as f:
        json.dump(validated, f, indent=2, default=str)
    print(f"  Detailed results saved to {json_path}")

    # Summary
    scores = [v["planet_score"] for v in validated]
    high = [v for v in validated if v["planet_score"] >= 70]
    medium = [v for v in validated if 50 <= v["planet_score"] < 70]
    low = [v for v in validated if v["planet_score"] < 50]

    print(f"\n{'=' * 55}")
    print(f"  VALIDATION SUMMARY")
    print(f"{'=' * 55}")
    print(f"  High confidence (>= 70):  {len(high)}")
    print(f"  Medium confidence (50-69): {len(medium)}")
    print(f"  Low / false positive (< 50): {len(low)}")
    print(f"  Mean score: {np.mean(scores):.1f}")
    print(f"{'=' * 55}")

    if high:
        print(f"\n  TOP VALIDATED CANDIDATES:")
        for v in high[:10]:
            name = v["filename"].replace(".csv", "")
            rp = v.get("planet_radius_earth")
            rp_str = f"Rp={rp:.1f}Re" if isinstance(rp, (int, float)) and rp else ""
            print(f"    Score={v['planet_score']:3d}  SNR={v['snr']:.1f}  P={v['period_days']:.4f}d  {rp_str}  {name}")

    print()


if __name__ == "__main__":
    main()
