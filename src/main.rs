use anyhow::{Context, Result};
use clap::Parser;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

// ============================================================================
// CLI
// ============================================================================

#[derive(Parser)]
#[command(
    name = "hunt",
    about = "🔭 Exoplanet Hunter — BLS Transit Detection in Rust"
)]
struct Cli {
    /// Directory containing light curve CSV files (time, flux, flux_err)
    #[arg(short, long)]
    input: PathBuf,

    /// Output JSON file for candidates
    #[arg(short, long, default_value = "candidates.json")]
    output: PathBuf,

    /// Minimum trial period in days
    #[arg(long, default_value_t = 0.5)]
    min_period: f64,

    /// Maximum trial period in days
    #[arg(long, default_value_t = 20.0)]
    max_period: f64,

    /// Number of trial periods to test
    #[arg(long, default_value_t = 10_000)]
    n_periods: usize,

    /// Number of phase bins
    #[arg(long, default_value_t = 200)]
    n_bins: usize,

    /// Minimum transit duration as fraction of period
    #[arg(long, default_value_t = 0.01)]
    min_duration_frac: f64,

    /// Maximum transit duration as fraction of period
    #[arg(long, default_value_t = 0.05)]
    max_duration_frac: f64,

    /// Minimum SNR to report as candidate
    #[arg(long, default_value_t = 6.0)]
    snr_threshold: f64,

    /// Number of threads (0 = auto)
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

// ============================================================================
// Data Structures
// ============================================================================

#[derive(Debug, Clone)]
struct LightCurve {
    filename: String,
    time: Vec<f64>,
    flux: Vec<f64>,
    flux_err: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct TransitCandidate {
    filename: String,
    period_days: f64,
    epoch: f64, // time of first transit
    duration_hours: f64,
    depth_ppm: f64, // transit depth in parts per million
    snr: f64,
    n_transits: usize,
    bls_power: f64,
    /// Estimated planet radius relative to star (Rp/Rs)
    radius_ratio: f64,
}

#[derive(Debug, Clone, Serialize)]
struct HuntReport {
    total_lightcurves: usize,
    candidates_found: usize,
    snr_threshold: f64,
    period_range: [f64; 2],
    candidates: Vec<TransitCandidate>,
}

// ============================================================================
// BLS Algorithm
// ============================================================================

/// Box-fitting Least Squares (Kovács, Zucker & Mazeh 2002)
///
/// For each trial period P:
///   1. Phase-fold the light curve
///   2. Bin the folded data
///   3. Slide a "box" (transit model) across all phases
///   4. For each box position and width, compute the BLS statistic:
///      SR = sqrt(n) * |signal_in_box - signal_out_box| / scatter
///   5. Record the best (highest SR) configuration
///
fn bls_search(
    lc: &LightCurve,
    periods: &[f64],
    n_bins: usize,
    min_dur_frac: f64,
    max_dur_frac: f64,
) -> Option<(f64, f64, f64, f64, f64)> {
    // (best_period, best_phase, best_duration_frac, best_depth, best_power)

    let n = lc.time.len();
    if n < 50 {
        return None; // not enough data
    }

    // Normalize flux: subtract median, divide by MAD
    let mut sorted_flux = lc.flux.clone();
    sorted_flux.sort_by_key(|&f| OrderedFloat(f));
    let median = sorted_flux[n / 2];

    let mut deviations: Vec<f64> = lc.flux.iter().map(|&f| (f - median).abs()).collect();
    deviations.sort_by_key(|&d| OrderedFloat(d));
    let mad = deviations[n / 2] * 1.4826; // scale to std

    if mad < 1e-10 {
        return None; // constant flux, no signal possible
    }

    let norm_flux: Vec<f64> = lc.flux.iter().map(|&f| (f - median) / mad).collect();

    // Weights from flux errors
    let weights: Vec<f64> = lc
        .flux_err
        .iter()
        .map(|&e| {
            let w = if e > 1e-10 { 1.0 / (e * e) } else { 1.0 };
            w
        })
        .collect();
    let total_weight: f64 = weights.iter().sum();

    let time_span = lc.time[n - 1] - lc.time[0];

    let mut best_power = f64::NEG_INFINITY;
    let mut best_period = 0.0;
    let mut best_phase = 0.0;
    let mut best_dur_frac = 0.0;
    let mut best_depth = 0.0;

    for &period in periods {
        if period > time_span * 0.5 {
            continue; // need at least 2 transits
        }

        // Phase-fold and bin
        let mut bin_sum = vec![0.0f64; n_bins];
        let mut bin_weight = vec![0.0f64; n_bins];
        let mut bin_count = vec![0usize; n_bins];

        for i in 0..n {
            let phase = ((lc.time[i] - lc.time[0]) % period) / period;
            let phase = if phase < 0.0 { phase + 1.0 } else { phase };
            let bin = (phase * n_bins as f64) as usize;
            let bin = bin.min(n_bins - 1);
            bin_sum[bin] += norm_flux[i] * weights[i];
            bin_weight[bin] += weights[i];
            bin_count[bin] += 1;
        }

        // Compute bin averages
        let _bin_avg: Vec<f64> = (0..n_bins)
            .map(|b| {
                if bin_weight[b] > 0.0 {
                    bin_sum[b] / bin_weight[b]
                } else {
                    0.0
                }
            })
            .collect();

        // Slide the box: try different transit start phases and durations
        let min_box = (min_dur_frac * n_bins as f64).max(1.0) as usize;
        let max_box = (max_dur_frac * n_bins as f64).max(min_box as f64 + 1.0) as usize;

        for box_width in min_box..=max_box.min(n_bins / 2) {
            for start in 0..n_bins {
                // Sum inside box
                let mut s_in = 0.0;
                let mut w_in = 0.0;
                let mut n_in = 0usize;

                for j in 0..box_width {
                    let bin = (start + j) % n_bins;
                    s_in += bin_sum[bin];
                    w_in += bin_weight[bin];
                    n_in += bin_count[bin];
                }

                if w_in < 1e-10 || n_in < 3 {
                    continue;
                }

                let w_out = total_weight - w_in;
                if w_out < 1e-10 {
                    continue;
                }

                let avg_in = s_in / w_in;
                let avg_out = (bin_sum.iter().sum::<f64>() - s_in) / w_out;
                let depth = avg_out - avg_in; // positive = dip

                // BLS power: SR statistic
                let power = depth * (w_in * w_out / total_weight).sqrt();

                if power > best_power {
                    best_power = power;
                    best_period = period;
                    best_phase = start as f64 / n_bins as f64;
                    best_dur_frac = box_width as f64 / n_bins as f64;
                    best_depth = depth;
                }
            }
        }
    }

    if best_power > f64::NEG_INFINITY {
        Some((
            best_period,
            best_phase,
            best_dur_frac,
            best_depth,
            best_power,
        ))
    } else {
        None
    }
}

/// Estimate SNR of a BLS detection by comparing to noise floor
fn estimate_snr(lc: &LightCurve, period: f64, phase: f64, dur_frac: f64) -> (f64, usize) {
    let n = lc.time.len();
    let t0 = lc.time[0];

    let mut in_transit_flux = Vec::new();
    let mut out_transit_flux = Vec::new();

    for i in 0..n {
        let p = ((lc.time[i] - t0) % period) / period;
        let p = if p < 0.0 { p + 1.0 } else { p };

        let in_box = if phase + dur_frac <= 1.0 {
            p >= phase && p < phase + dur_frac
        } else {
            p >= phase || p < (phase + dur_frac - 1.0)
        };

        if in_box {
            in_transit_flux.push(lc.flux[i]);
        } else {
            out_transit_flux.push(lc.flux[i]);
        }
    }

    if in_transit_flux.is_empty() || out_transit_flux.is_empty() {
        return (0.0, 0);
    }

    let mean_in: f64 = in_transit_flux.iter().sum::<f64>() / in_transit_flux.len() as f64;
    let mean_out: f64 = out_transit_flux.iter().sum::<f64>() / out_transit_flux.len() as f64;
    let depth = mean_out - mean_in;

    // Scatter of out-of-transit data
    let variance: f64 = out_transit_flux
        .iter()
        .map(|&f| (f - mean_out).powi(2))
        .sum::<f64>()
        / (out_transit_flux.len() as f64 - 1.0);
    let scatter = variance.sqrt();

    let n_in_transit = in_transit_flux.len();
    let snr = if scatter > 1e-10 {
        depth / (scatter / (n_in_transit as f64).sqrt())
    } else {
        0.0
    };

    // Count number of transits
    let time_span = lc.time[n - 1] - t0;
    let n_transits = (time_span / period).floor() as usize;

    (snr, n_transits.max(1))
}

// ============================================================================
// I/O
// ============================================================================

fn load_lightcurve(path: &Path) -> Result<LightCurve> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path)?;

    let mut time = Vec::new();
    let mut flux = Vec::new();
    let mut flux_err = Vec::new();

    for result in rdr.records() {
        let record = result?;
        if record.len() < 2 {
            continue;
        }

        let t: f64 = match record.get(0).unwrap_or("").parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let f: f64 = match record.get(1).unwrap_or("").parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let e: f64 = record.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.001);

        // Skip NaN or infinite values
        if t.is_finite() && f.is_finite() && e.is_finite() {
            time.push(t);
            flux.push(f);
            flux_err.push(e);
        }
    }

    let filename = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    Ok(LightCurve {
        filename,
        time,
        flux,
        flux_err,
    })
}

fn find_csv_files(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(dir).context("Cannot read input directory")? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(false, |e| e == "csv") {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

// ============================================================================
// Main
// ============================================================================

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    println!("\n🔭 Exoplanet Hunter v0.1.0");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!(
        "  Period range:  {:.2} – {:.2} days",
        cli.min_period, cli.max_period
    );
    println!("  Trial periods: {}", cli.n_periods);
    println!("  Phase bins:    {}", cli.n_bins);
    println!("  SNR threshold: {:.1}", cli.snr_threshold);
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Generate trial periods (log-spaced for even coverage)
    let log_min = cli.min_period.ln();
    let log_max = cli.max_period.ln();
    let periods: Vec<f64> = (0..cli.n_periods)
        .map(|i| {
            let frac = i as f64 / (cli.n_periods - 1) as f64;
            (log_min + frac * (log_max - log_min)).exp()
        })
        .collect();

    // Find all light curve files
    let files = find_csv_files(&cli.input)?;
    println!("📁 Found {} light curve files\n", files.len());

    if files.is_empty() {
        println!("⚠️  No CSV files found in {:?}", cli.input);
        println!("   Expected format: time,flux,flux_err");
        return Ok(());
    }

    // Progress bar
    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("  Hunting {bar:40.cyan/blue} {pos}/{len} stars [{elapsed_precise}] {msg}")
            .unwrap()
            .progress_chars("█▉▊▋▌▍▎▏ "),
    );

    // Process all light curves in parallel
    let counter = AtomicU64::new(0);
    let candidates: Vec<TransitCandidate> = files
        .par_iter()
        .filter_map(|path| {
            let done = counter.fetch_add(1, Ordering::Relaxed);
            pb.set_position(done + 1);
            let lc = match load_lightcurve(path) {
                Ok(lc) => lc,
                Err(_) => return None,
            };

            if lc.time.len() < 100 {
                return None; // too few points
            }

            // Run BLS
            let result = bls_search(
                &lc,
                &periods,
                cli.n_bins,
                cli.min_duration_frac,
                cli.max_duration_frac,
            )?;

            let (period, phase, dur_frac, depth, power) = result;

            // Compute SNR
            let (snr, n_transits) = estimate_snr(&lc, period, phase, dur_frac);

            if snr >= cli.snr_threshold && n_transits >= 2 && depth > 0.0 {
                let duration_hours = dur_frac * period * 24.0;
                let depth_ppm = depth * 1e6; // assuming normalized flux
                let radius_ratio = depth.abs().sqrt(); // Rp/Rs ≈ sqrt(depth)
                let epoch = lc.time[0] + phase * period;

                Some(TransitCandidate {
                    filename: lc.filename,
                    period_days: period,
                    epoch,
                    duration_hours,
                    depth_ppm,
                    snr,
                    n_transits,
                    bls_power: power,
                    radius_ratio,
                })
            } else {
                None
            }
        })
        .collect();

    pb.finish_with_message("done!");

    // Sort by SNR descending
    let mut candidates = candidates;
    candidates.sort_by(|a, b| {
        b.snr
            .partial_cmp(&a.snr)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Print results
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

    // Write JSON report
    let report = HuntReport {
        total_lightcurves: files.len(),
        candidates_found: candidates.len(),
        snr_threshold: cli.snr_threshold,
        period_range: [cli.min_period, cli.max_period],
        candidates,
    };

    let json = serde_json::to_string_pretty(&report)?;
    fs::write(&cli.output, &json)?;
    println!("\n📄 Report saved to {:?}\n", cli.output);

    Ok(())
}
