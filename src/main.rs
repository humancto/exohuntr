use anyhow::Result;
use clap::{Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::Serialize;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use exoplanet_hunter::{bls, crossmatch, io, validate};

// ============================================================================
// CLI
// ============================================================================

#[derive(Parser)]
#[command(
    name = "hunt",
    about = "Exoplanet Hunter — BLS Transit Detection & Validation in Rust"
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    // Top-level flags for backward compatibility (search mode)
    /// Directory containing light curve CSV files
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output JSON file for candidates
    #[arg(short, long, default_value = "candidates.json")]
    output: PathBuf,

    #[arg(long, default_value_t = 0.5)]
    min_period: f64,

    #[arg(long, default_value_t = 20.0)]
    max_period: f64,

    #[arg(long, default_value_t = 10_000)]
    n_periods: usize,

    #[arg(long, default_value_t = 200)]
    n_bins: usize,

    #[arg(long, default_value_t = 0.01)]
    min_duration_frac: f64,

    #[arg(long, default_value_t = 0.05)]
    max_duration_frac: f64,

    #[arg(long, default_value_t = 6.0)]
    snr_threshold: f64,

    /// Number of threads (0 = auto)
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

#[derive(Subcommand)]
enum Commands {
    /// Run BLS transit search on light curves
    Search {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value = "candidates.json")]
        output: PathBuf,
        #[arg(long, default_value_t = 0.5)]
        min_period: f64,
        #[arg(long, default_value_t = 20.0)]
        max_period: f64,
        #[arg(long, default_value_t = 10_000)]
        n_periods: usize,
        #[arg(long, default_value_t = 200)]
        n_bins: usize,
        #[arg(long, default_value_t = 0.01)]
        min_duration_frac: f64,
        #[arg(long, default_value_t = 0.05)]
        max_duration_frac: f64,
        #[arg(long, default_value_t = 6.0)]
        snr_threshold: f64,
        #[arg(long, default_value_t = 0)]
        threads: usize,
    },
    /// Validate BLS candidates with false-positive tests
    Validate {
        /// Input candidates JSON from search
        #[arg(short, long, default_value = "candidates.json")]
        input: PathBuf,
        /// Directory containing light curve CSV files
        #[arg(short, long, default_value = "data/lightcurves")]
        lightcurves: PathBuf,
        /// Output directory for validation results
        #[arg(short, long, default_value = "results")]
        output: PathBuf,
        /// Number of threads (0 = auto)
        #[arg(long, default_value_t = 0)]
        threads: usize,
    },
    /// Cross-match candidates against known exoplanet catalog
    Crossmatch {
        /// Input candidates JSON from search
        #[arg(short, long, default_value = "candidates.json")]
        input: PathBuf,
        /// Path to confirmed exoplanet catalog CSV
        #[arg(short, long, default_value = "data/lightcurves/confirmed_exoplanets.csv")]
        catalog: PathBuf,
        /// Output CSV for cross-match results
        #[arg(short, long, default_value = "results/crossmatch_results.csv")]
        output: PathBuf,
    },
}

// ============================================================================
// Data Structures
// ============================================================================

#[derive(Debug, Clone, Serialize, serde::Deserialize)]
struct TransitCandidate {
    filename: String,
    period_days: f64,
    epoch: f64,
    duration_hours: f64,
    depth_ppm: f64,
    snr: f64,
    n_transits: usize,
    bls_power: f64,
    radius_ratio: f64,
}

#[derive(Debug, Clone, Serialize, serde::Deserialize)]
struct HuntReport {
    total_lightcurves: usize,
    candidates_found: usize,
    snr_threshold: f64,
    period_range: [f64; 2],
    candidates: Vec<TransitCandidate>,
}

// ============================================================================
// Search
// ============================================================================

fn run_search(
    input: &PathBuf,
    output: &PathBuf,
    min_period: f64,
    max_period: f64,
    n_periods: usize,
    n_bins: usize,
    min_duration_frac: f64,
    max_duration_frac: f64,
    snr_threshold: f64,
    threads: usize,
) -> Result<()> {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }

    println!("\n🔭 Exoplanet Hunter v0.2.0");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  Period range:  {:.2} – {:.2} days", min_period, max_period);
    println!("  Trial periods: {}", n_periods);
    println!("  Phase bins:    {}", n_bins);
    println!("  SNR threshold: {:.1}", snr_threshold);
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    let periods = bls::generate_periods(min_period, max_period, n_periods);
    let files = io::find_csv_files(input)?;
    println!("📁 Found {} light curve files\n", files.len());

    if files.is_empty() {
        println!("⚠️  No CSV files found in {:?}", input);
        return Ok(());
    }

    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("  Hunting {bar:40.cyan/blue} {pos}/{len} stars [{elapsed_precise}] {msg}")
            .unwrap()
            .progress_chars("█▉▊▋▌▍▎▏ "),
    );

    let counter = AtomicU64::new(0);
    let candidates: Vec<TransitCandidate> = files
        .par_iter()
        .filter_map(|path| {
            let done = counter.fetch_add(1, Ordering::Relaxed);
            pb.set_position(done + 1);
            let lc = io::load_lightcurve(path).ok()?;

            if lc.time.len() < 100 {
                return None;
            }

            let result = bls::bls_search(&lc, &periods, n_bins, min_duration_frac, max_duration_frac)?;

            let (snr, n_transits) = bls::estimate_snr(&lc, result.period, result.phase, result.duration_frac);

            if snr >= snr_threshold && n_transits >= 2 && result.depth > 0.0 {
                let duration_hours = result.duration_frac * result.period * 24.0;
                let depth_ppm = result.depth * 1e6;
                let radius_ratio = result.depth.abs().sqrt();
                let epoch = lc.time[0] + result.phase * result.period;

                Some(TransitCandidate {
                    filename: lc.filename,
                    period_days: result.period,
                    epoch,
                    duration_hours,
                    depth_ppm,
                    snr,
                    n_transits,
                    bls_power: result.power,
                    radius_ratio,
                })
            } else {
                None
            }
        })
        .collect();

    pb.finish_with_message("done!");

    let mut candidates = candidates;
    candidates.sort_by(|a, b| {
        b.snr
            .partial_cmp(&a.snr)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    println!("\n\n🎯 CANDIDATES FOUND: {}", candidates.len());
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    if candidates.is_empty() {
        println!("  No candidates above SNR threshold. Try lowering --snr-threshold.");
    } else {
        println!(
            "  {:>4} {:>30} {:>10} {:>8} {:>10} {:>8} {:>8}",
            "Rank", "File", "Period(d)", "SNR", "Depth(ppm)", "Transits", "Rp/Rs"
        );
        println!("  {}", "─".repeat(88));

        for (i, c) in candidates.iter().enumerate().take(50) {
            let short_name: String = c.filename.chars().take(28).collect();
            println!(
                "  {:>4} {:>30} {:>10.4} {:>8.1} {:>10.0} {:>8} {:>8.4}",
                i + 1,
                short_name,
                c.period_days,
                c.snr,
                c.depth_ppm,
                c.n_transits,
                c.radius_ratio,
            );
        }
    }

    let report = HuntReport {
        total_lightcurves: files.len(),
        candidates_found: candidates.len(),
        snr_threshold,
        period_range: [min_period, max_period],
        candidates,
    };

    let json = serde_json::to_string_pretty(&report)?;
    fs::write(output, &json)?;
    println!("\n📄 Report saved to {:?}\n", output);

    Ok(())
}

// ============================================================================
// Validate
// ============================================================================

fn run_validate(
    input: &PathBuf,
    lightcurves: &PathBuf,
    output: &PathBuf,
    threads: usize,
) -> Result<()> {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }

    println!("\nExohuntr — Rust Validation Pipeline");
    println!("{}", "=".repeat(55));

    let data = fs::read_to_string(input)?;
    let report: HuntReport = serde_json::from_str(&data)?;
    println!("  Loaded {} candidates\n", report.candidates.len());

    let results: Vec<validate::ValidationResult> = report
        .candidates
        .par_iter()
        .filter_map(|c| {
            let filepath = lightcurves.join(&c.filename);
            let lc = io::load_lightcurve(&filepath).ok()?;

            Some(validate::validate_candidate(
                &lc,
                c.period_days,
                c.epoch,
                c.duration_hours,
                c.snr,
                c.depth_ppm,
                c.radius_ratio,
                c.n_transits,
                c.bls_power,
                None, // reference period from ExoFOP handled in Python
            ))
        })
        .collect();

    // Output
    fs::create_dir_all(output)?;
    let json_path = output.join("validation_results.json");
    let json = serde_json::to_string_pretty(&results)?;
    fs::write(&json_path, &json)?;

    // Summary
    let high: Vec<_> = results.iter().filter(|r| r.planet_score >= 70).collect();
    let medium: Vec<_> = results
        .iter()
        .filter(|r| r.planet_score >= 50 && r.planet_score < 70)
        .collect();
    let low: Vec<_> = results.iter().filter(|r| r.planet_score < 50).collect();

    println!("{}", "=".repeat(55));
    println!("  VALIDATION SUMMARY");
    println!("{}", "=".repeat(55));
    println!("  Total validated:           {}", results.len());
    println!("  High confidence (>= 70):   {}", high.len());
    println!("  Medium confidence (50-69):  {}", medium.len());
    println!("  Low / false positive (< 50): {}", low.len());
    println!("{}", "=".repeat(55));

    if !high.is_empty() {
        println!("\n  TOP VALIDATED CANDIDATES:");
        for v in high.iter().take(10) {
            println!(
                "    Score={:3}  SNR={:.1}  P={:.4}d  Rp/Rs={:.4}  {}",
                v.planet_score, v.snr, v.period_days, v.radius_ratio, v.filename
            );
        }
    }

    println!("\n📄 Results saved to {:?}\n", json_path);
    Ok(())
}

// ============================================================================
// Crossmatch
// ============================================================================

fn run_crossmatch(input: &PathBuf, catalog: &PathBuf, output: &PathBuf) -> Result<()> {
    println!("\nExohuntr — Cross-Match Pipeline");
    println!("{}", "=".repeat(55));

    let data = fs::read_to_string(input)?;
    let report: HuntReport = serde_json::from_str(&data)?;
    println!("  Loaded {} candidates", report.candidates.len());

    let entries = crossmatch::load_catalog(catalog)?;
    println!("  Loaded {} catalog entries", entries.len());

    let index = crossmatch::CatalogIndex::new(entries);

    let candidate_tuples: Vec<(String, f64, f64)> = report
        .candidates
        .iter()
        .map(|c| (c.filename.clone(), c.period_days, c.snr))
        .collect();

    let results = crossmatch::crossmatch_candidates(&candidate_tuples, &index);

    let known_count = results.iter().filter(|r| r.status == "KNOWN").count();
    let new_count = results.iter().filter(|r| r.status == "POTENTIALLY NEW").count();

    // Write CSV
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut wtr = csv::Writer::from_path(output)?;
    for r in &results {
        wtr.serialize(r)?;
    }
    wtr.flush()?;

    println!("  Known: {} | Potentially new: {}", known_count, new_count);
    println!("📄 Results saved to {:?}\n", output);
    Ok(())
}

// ============================================================================
// Main
// ============================================================================

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Search {
            input,
            output,
            min_period,
            max_period,
            n_periods,
            n_bins,
            min_duration_frac,
            max_duration_frac,
            snr_threshold,
            threads,
        }) => run_search(
            &input,
            &output,
            min_period,
            max_period,
            n_periods,
            n_bins,
            min_duration_frac,
            max_duration_frac,
            snr_threshold,
            threads,
        ),
        Some(Commands::Validate {
            input,
            lightcurves,
            output,
            threads,
        }) => run_validate(&input, &lightcurves, &output, threads),
        Some(Commands::Crossmatch {
            input,
            catalog,
            output,
        }) => run_crossmatch(&input, &catalog, &output),
        None => {
            // Backward compatibility: if --input is provided, run search
            if let Some(input) = cli.input {
                run_search(
                    &input,
                    &cli.output,
                    cli.min_period,
                    cli.max_period,
                    cli.n_periods,
                    cli.n_bins,
                    cli.min_duration_frac,
                    cli.max_duration_frac,
                    cli.snr_threshold,
                    cli.threads,
                )
            } else {
                println!("Usage: hunt <command> or hunt -i <input_dir>");
                println!("Commands: search, validate, crossmatch");
                println!("Run 'hunt --help' for details.");
                Ok(())
            }
        }
    }
}
