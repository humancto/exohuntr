#!/usr/bin/env python3
"""
📊 analyze_candidates.py — Analyze and visualize BLS transit detection results.

Takes the candidates.json output from the Rust BLS engine and generates:
  1. Phase-folded light curve plots for each candidate
  2. Period vs depth scatter plot
  3. Cross-reference with known exoplanet catalogs
  4. A final report ranking candidates by discovery potential

Usage:
    python analyze_candidates.py --input candidates.json --lightcurves data/lightcurves/
    python analyze_candidates.py --input candidates.json --crossmatch  # check if any are NEW
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================================
# Plot styling
# ============================================================================

COLORS = {
    "bg": "#0a0e1a",
    "fg": "#e2e8f0",
    "accent": "#22d3ee",
    "orange": "#fb923c",
    "purple": "#a78bfa",
    "green": "#34d399",
    "grid": "#1e293b",
    "transit": "#ef4444",
}

def setup_style():
    plt.rcParams.update({
        "figure.facecolor": COLORS["bg"],
        "axes.facecolor": COLORS["bg"],
        "axes.edgecolor": COLORS["grid"],
        "axes.labelcolor": COLORS["fg"],
        "text.color": COLORS["fg"],
        "xtick.color": COLORS["fg"],
        "ytick.color": COLORS["fg"],
        "grid.color": COLORS["grid"],
        "grid.alpha": 0.3,
        "font.family": "sans-serif",
        "font.size": 11,
    })


# ============================================================================
# Phase-folded light curve plot
# ============================================================================

def plot_phase_folded(candidate: dict, lc_dir: Path, output_dir: Path):
    """Create a phase-folded light curve plot for a candidate."""
    filepath = lc_dir / candidate["filename"]
    if not filepath.exists():
        return None

    df = pd.read_csv(filepath)
    time = df["time"].values
    flux = df["flux"].values

    period = candidate["period_days"]
    epoch = candidate["epoch"]

    # Phase fold
    phase = ((time - epoch) / period) % 1.0
    phase[phase > 0.5] -= 1.0  # center transit at 0

    # Bin the folded data
    n_bins = 100
    bin_edges = np.linspace(-0.5, 0.5, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_flux = np.zeros(n_bins)
    bin_count = np.zeros(n_bins)

    for i in range(len(phase)):
        b = int((phase[i] + 0.5) * n_bins)
        b = min(b, n_bins - 1)
        bin_flux[b] += flux[i]
        bin_count[b] += 1

    mask = bin_count > 0
    bin_flux[mask] /= bin_count[mask]

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), height_ratios=[3, 1])
    fig.suptitle(
        f"🔭 {candidate['filename']}  |  P = {period:.4f} d  |  SNR = {candidate['snr']:.1f}",
        fontsize=13, fontweight="bold", y=0.98,
    )

    # Top: phase-folded with individual points
    ax1.scatter(phase, flux, s=1, alpha=0.15, color=COLORS["fg"], rasterized=True)
    ax1.plot(bin_centers[mask], bin_flux[mask], "o-", color=COLORS["accent"],
             markersize=4, linewidth=1.5, label="Binned", zorder=5)

    # Mark transit region
    dur_phase = candidate["duration_hours"] / (period * 24)
    ax1.axvspan(-dur_phase / 2, dur_phase / 2, alpha=0.1, color=COLORS["transit"], label="Transit window")
    ax1.axvline(0, color=COLORS["transit"], alpha=0.3, linestyle="--", linewidth=0.8)

    ax1.set_ylabel("Normalized Flux")
    ax1.set_xlim(-0.5, 0.5)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.15)

    # Bottom: zoom on transit
    zoom_range = dur_phase * 3
    zoom_mask = (phase > -zoom_range) & (phase < zoom_range)
    ax2.scatter(phase[zoom_mask], flux[zoom_mask], s=3, alpha=0.3, color=COLORS["fg"], rasterized=True)
    binmask = mask & (bin_centers > -zoom_range) & (bin_centers < zoom_range)
    ax2.plot(bin_centers[binmask], bin_flux[binmask], "o-", color=COLORS["orange"],
             markersize=5, linewidth=2, zorder=5)
    ax2.axvspan(-dur_phase / 2, dur_phase / 2, alpha=0.15, color=COLORS["transit"])
    ax2.set_xlabel("Orbital Phase")
    ax2.set_ylabel("Flux (zoom)")
    ax2.set_xlim(-zoom_range, zoom_range)
    ax2.grid(True, alpha=0.15)

    # Annotation box
    info = (
        f"Depth: {candidate['depth_ppm']:.0f} ppm\n"
        f"Duration: {candidate['duration_hours']:.2f} h\n"
        f"Rp/Rs: {candidate['radius_ratio']:.4f}\n"
        f"Transits: {candidate['n_transits']}"
    )
    ax1.text(0.02, 0.02, info, transform=ax1.transAxes, fontsize=9,
             verticalalignment="bottom", fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor=COLORS["grid"], alpha=0.8))

    plt.tight_layout()
    safe_name = candidate["filename"].replace(".csv", "")
    outpath = output_dir / f"phase_fold_{safe_name}.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outpath


# ============================================================================
# Overview plots
# ============================================================================

def plot_candidate_overview(candidates: list, output_dir: Path):
    """Create overview scatter plots of all candidates."""
    if not candidates:
        return

    periods = [c["period_days"] for c in candidates]
    depths = [c["depth_ppm"] for c in candidates]
    snrs = [c["snr"] for c in candidates]
    radii = [c["radius_ratio"] for c in candidates]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"🔭 Exoplanet Hunter — {len(candidates)} Candidates", fontsize=14, fontweight="bold")

    # 1. Period vs Depth
    ax = axes[0]
    sc = ax.scatter(periods, depths, c=snrs, cmap="plasma", s=40, alpha=0.8, edgecolors="white", linewidths=0.3)
    ax.set_xlabel("Orbital Period (days)")
    ax.set_ylabel("Transit Depth (ppm)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.15)
    plt.colorbar(sc, ax=ax, label="SNR", shrink=0.8)

    # 2. Period vs Radius Ratio
    ax = axes[1]
    sc = ax.scatter(periods, radii, c=snrs, cmap="plasma", s=40, alpha=0.8, edgecolors="white", linewidths=0.3)
    ax.set_xlabel("Orbital Period (days)")
    ax.set_ylabel("Rp/Rs (Planet/Star radius ratio)")
    ax.set_xscale("log")
    ax.grid(True, alpha=0.15)

    # Reference lines for planet sizes (assuming Sun-like star)
    for label, ratio in [("Earth", 0.009), ("Neptune", 0.036), ("Jupiter", 0.1)]:
        ax.axhline(ratio, color=COLORS["green"], alpha=0.3, linestyle="--", linewidth=0.8)
        ax.text(periods[0] if periods else 1, ratio * 1.1, label, fontsize=8, color=COLORS["green"], alpha=0.6)

    plt.colorbar(sc, ax=ax, label="SNR", shrink=0.8)

    # 3. SNR histogram
    ax = axes[2]
    ax.hist(snrs, bins=30, color=COLORS["accent"], alpha=0.7, edgecolor=COLORS["bg"])
    ax.axvline(np.median(snrs), color=COLORS["orange"], linestyle="--", label=f"Median: {np.median(snrs):.1f}")
    ax.set_xlabel("Signal-to-Noise Ratio")
    ax.set_ylabel("Count")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.15)

    plt.tight_layout()
    fig.savefig(output_dir / "candidate_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   📊 Saved overview plot")


# ============================================================================
# Cross-matching
# ============================================================================

def crossmatch_known_planets(candidates: list, catalog_path: Path) -> pd.DataFrame:
    """Cross-match candidates with the known exoplanet catalog."""
    if not catalog_path.exists():
        print("   ⚠️  No catalog file found. Run downloader with --catalog first.")
        return pd.DataFrame()

    catalog = pd.read_csv(catalog_path)
    results = []

    for c in candidates:
        # Extract target identifier from filename
        fname = c["filename"]

        # Try to match by TIC ID or KIC ID
        matched = False
        for _, planet in catalog.iterrows():
            if str(planet.get("hostname", "")).replace(" ", "_") in fname:
                results.append({
                    "candidate_file": fname,
                    "candidate_period": c["period_days"],
                    "candidate_snr": c["snr"],
                    "known_planet": planet.get("pl_name", ""),
                    "known_period": planet.get("pl_orbper", None),
                    "known_radius_earth": planet.get("pl_rade", None),
                    "status": "KNOWN",
                })
                matched = True
                break

        if not matched:
            results.append({
                "candidate_file": fname,
                "candidate_period": c["period_days"],
                "candidate_snr": c["snr"],
                "known_planet": "",
                "known_period": None,
                "known_radius_earth": None,
                "status": "⭐ POTENTIALLY NEW",
            })

    df = pd.DataFrame(results)
    return df


# ============================================================================
# Report generation
# ============================================================================

def generate_report(report: dict, crossmatch_df: pd.DataFrame, output_dir: Path):
    """Generate a markdown report."""
    candidates = report["candidates"]

    lines = [
        "# 🔭 Exoplanet Hunter — Analysis Report\n",
        f"**Total light curves analyzed:** {report['total_lightcurves']}",
        f"**Candidates found:** {report['candidates_found']}",
        f"**SNR threshold:** {report['snr_threshold']}",
        f"**Period range:** {report['period_range'][0]:.2f} – {report['period_range'][1]:.2f} days\n",
        "---\n",
        "## Top Candidates\n",
        "| Rank | File | Period (d) | SNR | Depth (ppm) | Rp/Rs | Status |",
        "|------|------|-----------|-----|-------------|-------|--------|",
    ]

    for i, c in enumerate(candidates[:30]):
        status = "?"
        if not crossmatch_df.empty:
            match = crossmatch_df[crossmatch_df["candidate_file"] == c["filename"]]
            if len(match) > 0:
                status = match.iloc[0]["status"]

        lines.append(
            f"| {i+1} | `{c['filename'][:35]}` | {c['period_days']:.4f} | "
            f"{c['snr']:.1f} | {c['depth_ppm']:.0f} | {c['radius_ratio']:.4f} | {status} |"
        )

    # Potentially new discoveries
    if not crossmatch_df.empty:
        new_ones = crossmatch_df[crossmatch_df["status"] == "⭐ POTENTIALLY NEW"]
        if len(new_ones) > 0:
            lines.extend([
                "\n---\n",
                f"## ⭐ Potentially New Discoveries ({len(new_ones)})\n",
                "These candidates were NOT found in the confirmed exoplanet catalog.\n",
                "**⚠️ Important:** These need further validation before claiming discovery:",
                "- Check for eclipsing binaries (secondary eclipse, odd/even transit depth differences)",
                "- Look for centroid shifts (signal might be from a background star)",
                "- Cross-check with TESS TOI list and community follow-up observations",
                "- Submit to [ExoFOP](https://exofop.ipac.caltech.edu/tess/) for community vetting\n",
            ])

            for _, row in new_ones.iterrows():
                lines.append(f"- **{row['candidate_file']}** — P={row['candidate_period']:.4f}d, SNR={row['candidate_snr']:.1f}")

    lines.extend([
        "\n---\n",
        "## Phase-Folded Light Curves\n",
        "See the `plots/` directory for phase-folded light curve visualizations of each candidate.\n",
        "\n---\n",
        "## Next Steps\n",
        "1. **Visual inspection**: Review phase-folded plots for transit-like shapes",
        "2. **False positive tests**: Check for V-shaped eclipses (binaries), centroid motion",
        "3. **Period validation**: Look for secondary eclipse at phase 0.5",
        "4. **Cross-reference**: Check ExoFOP, SIMBAD, and recent literature",
        "5. **Follow-up**: If promising, submit to ExoFOP or contact a professional astronomer",
        "6. **Publish**: Write up findings, share on GitHub, submit to citizen science programs\n",
        "\n*Generated by Exoplanet Hunter (github.com/your-username/exoplanet-hunter)*",
    ])

    report_path = output_dir / "REPORT.md"
    report_path.write_text("\n".join(lines))
    print(f"   📝 Report saved to {report_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="📊 Analyze BLS transit detection results")
    parser.add_argument("--input", default="candidates.json", help="BLS results JSON")
    parser.add_argument("--lightcurves", default="data/lightcurves", help="Light curve directory")
    parser.add_argument("--output", default="results", help="Output directory for plots and report")
    parser.add_argument("--crossmatch", action="store_true", help="Cross-match with known exoplanet catalog")
    parser.add_argument("--catalog", default="data/lightcurves/confirmed_exoplanets.csv", help="Catalog CSV path")
    parser.add_argument("--top-n", type=int, default=20, help="Generate phase-fold plots for top N candidates")

    args = parser.parse_args()

    setup_style()
    output_dir = Path(args.output)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load results
    print("\n📊 Exoplanet Hunter — Analysis")
    print("━" * 50)

    with open(args.input) as f:
        report = json.load(f)

    candidates = report["candidates"]
    print(f"   Loaded {len(candidates)} candidates from {args.input}")

    if not candidates:
        print("   No candidates to analyze. Try lowering the SNR threshold.")
        sys.exit(0)

    # Overview plots
    plot_candidate_overview(candidates, plots_dir)

    # Phase-folded plots for top candidates
    lc_dir = Path(args.lightcurves)
    print(f"\n   Generating phase-fold plots for top {args.top_n} candidates...")
    for c in candidates[: args.top_n]:
        result = plot_phase_folded(c, lc_dir, plots_dir)
        if result:
            print(f"   ✅ {c['filename']} (SNR={c['snr']:.1f})")

    # Cross-match
    crossmatch_df = pd.DataFrame()
    if args.crossmatch:
        print("\n   🔍 Cross-matching with known exoplanet catalog...")
        crossmatch_df = crossmatch_known_planets(candidates, Path(args.catalog))
        if not crossmatch_df.empty:
            new_count = len(crossmatch_df[crossmatch_df["status"] == "⭐ POTENTIALLY NEW"])
            known_count = len(crossmatch_df[crossmatch_df["status"] == "KNOWN"])
            print(f"   Known: {known_count} | Potentially new: {new_count}")

            crossmatch_df.to_csv(output_dir / "crossmatch_results.csv", index=False)

    # Generate report
    generate_report(report, crossmatch_df, output_dir)

    print(f"\n{'━' * 50}")
    print(f"🏁 Analysis complete!")
    print(f"   📂 Results: {output_dir.resolve()}")
    print(f"   📊 Plots:   {plots_dir.resolve()}")
    print(f"   📝 Report:  {(output_dir / 'REPORT.md').resolve()}")


if __name__ == "__main__":
    main()
