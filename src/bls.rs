use crate::io::LightCurve;
use ordered_float::OrderedFloat;

/// Result of a BLS search on a single light curve.
#[derive(Debug, Clone)]
pub struct BlsResult {
    pub period: f64,
    pub phase: f64,
    pub duration_frac: f64,
    pub depth: f64,
    pub power: f64,
}

/// Box-fitting Least Squares (Kovacs, Zucker & Mazeh 2002).
///
/// For each trial period, phase-folds the light curve into bins, then slides
/// a box-shaped transit model to find the configuration with maximum BLS power.
///
/// Returns `None` if the light curve has fewer than 50 points or constant flux.
pub fn bls_search(
    lc: &LightCurve,
    periods: &[f64],
    n_bins: usize,
    min_dur_frac: f64,
    max_dur_frac: f64,
) -> Option<BlsResult> {
    let n = lc.time.len();
    if n < 50 {
        return None;
    }

    // Normalize flux: subtract median, divide by MAD
    let mut sorted_flux = lc.flux.clone();
    sorted_flux.sort_by_key(|&f| OrderedFloat(f));
    let median = sorted_flux[n / 2];

    let mut deviations: Vec<f64> = lc.flux.iter().map(|&f| (f - median).abs()).collect();
    deviations.sort_by_key(|&d| OrderedFloat(d));
    let mad = deviations[n / 2] * 1.4826;

    if mad < 1e-10 {
        return None;
    }

    let norm_flux: Vec<f64> = lc.flux.iter().map(|&f| (f - median) / mad).collect();

    let weights: Vec<f64> = lc
        .flux_err
        .iter()
        .map(|&e| if e > 1e-10 { 1.0 / (e * e) } else { 1.0 })
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
            continue;
        }

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

        let min_box = (min_dur_frac * n_bins as f64).max(1.0) as usize;
        let max_box = (max_dur_frac * n_bins as f64).max(min_box as f64 + 1.0) as usize;

        let total_sum: f64 = bin_sum.iter().sum();

        for box_width in min_box..=max_box.min(n_bins / 2) {
            for start in 0..n_bins {
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
                let avg_out = (total_sum - s_in) / w_out;
                let depth = avg_out - avg_in;

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
        Some(BlsResult {
            period: best_period,
            phase: best_phase,
            duration_frac: best_dur_frac,
            depth: best_depth,
            power: best_power,
        })
    } else {
        None
    }
}

/// Estimate SNR of a BLS detection by comparing in-transit to out-of-transit flux.
///
/// Returns (snr, n_transits). The SNR is computed as the transit depth divided
/// by the per-point uncertainty of the in-transit mean.
pub fn estimate_snr(lc: &LightCurve, period: f64, phase: f64, dur_frac: f64) -> (f64, usize) {
    let n = lc.time.len();
    if n == 0 {
        return (0.0, 0);
    }
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

    if in_transit_flux.is_empty() || out_transit_flux.len() < 2 {
        return (0.0, 0);
    }

    let mean_in: f64 = in_transit_flux.iter().sum::<f64>() / in_transit_flux.len() as f64;
    let mean_out: f64 = out_transit_flux.iter().sum::<f64>() / out_transit_flux.len() as f64;
    let depth = mean_out - mean_in;

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

    let time_span = if n > 1 { lc.time[n - 1] - t0 } else { 0.0 };
    let n_transits = if period > 0.0 {
        (time_span / period).floor() as usize
    } else {
        0
    };

    (snr, n_transits.max(1))
}

/// Generate log-spaced trial periods between min_period and max_period.
///
/// # Examples
///
/// ```
/// use exoplanet_hunter::bls::generate_periods;
///
/// let periods = generate_periods(1.0, 10.0, 100);
/// assert_eq!(periods.len(), 100);
/// assert!((periods[0] - 1.0).abs() < 1e-10);
/// ```
pub fn generate_periods(min_period: f64, max_period: f64, n_periods: usize) -> Vec<f64> {
    let log_min = min_period.ln();
    let log_max = max_period.ln();
    (0..n_periods)
        .map(|i| {
            let frac = i as f64 / (n_periods - 1) as f64;
            (log_min + frac * (log_max - log_min)).exp()
        })
        .collect()
}

/// Compute phase array for a light curve given period and epoch.
///
/// # Examples
///
/// ```
/// use exoplanet_hunter::bls::compute_phases;
///
/// let phases = compute_phases(&[0.0, 1.0, 2.0], 2.0, 0.0);
/// assert!((phases[0] - 0.0).abs() < 1e-10);
/// assert!((phases[1] - 0.5).abs() < 1e-10);
/// ```
pub fn compute_phases(time: &[f64], period: f64, epoch: f64) -> Vec<f64> {
    time.iter()
        .map(|&t| {
            let p = ((t - epoch) / period) % 1.0;
            if p < 0.0 { p + 1.0 } else { p }
        })
        .collect()
}

/// Compute median of a slice. Returns 0.0 for empty input.
///
/// # Examples
///
/// ```
/// use exoplanet_hunter::bls::median;
///
/// assert!((median(&[3.0, 1.0, 2.0]) - 2.0).abs() < 1e-10);
/// assert!((median(&[4.0, 1.0, 3.0, 2.0]) - 2.5).abs() < 1e-10);
/// assert_eq!(median(&[]), 0.0);
/// ```
pub fn median(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by_key(|&v| OrderedFloat(v));
    let n = sorted.len();
    if n.is_multiple_of(2) {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Create a synthetic light curve with an injected box transit and small noise.
    ///
    /// Adds deterministic pseudo-noise so MAD is nonzero (required for BLS normalization).
    fn make_transit_lc(
        n_points: usize,
        period: f64,
        phase_center: f64,
        duration_frac: f64,
        depth: f64,
        time_span: f64,
    ) -> LightCurve {
        let mut time = Vec::with_capacity(n_points);
        let mut flux = Vec::with_capacity(n_points);
        let mut flux_err = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let t = i as f64 / n_points as f64 * time_span;
            let p = ((t % period) / period) % 1.0;

            let half_dur = duration_frac / 2.0;
            let in_transit = if phase_center - half_dur < 0.0 {
                p < phase_center + half_dur || p > 1.0 + phase_center - half_dur
            } else if phase_center + half_dur > 1.0 {
                p > phase_center - half_dur || p < phase_center + half_dur - 1.0
            } else {
                p >= phase_center - half_dur && p <= phase_center + half_dur
            };

            // Deterministic pseudo-noise: small sinusoidal variation
            let noise = 0.0001 * ((i as f64 * 0.1).sin() + (i as f64 * 0.037).cos());
            let f = if in_transit {
                1.0 - depth + noise
            } else {
                1.0 + noise
            };
            time.push(t);
            flux.push(f);
            flux_err.push(0.001);
        }

        LightCurve {
            filename: "synthetic.csv".to_string(),
            time,
            flux,
            flux_err,
        }
    }

    /// Create a constant-flux light curve (no transit).
    fn make_flat_lc(n_points: usize, time_span: f64) -> LightCurve {
        let time: Vec<f64> = (0..n_points)
            .map(|i| i as f64 / n_points as f64 * time_span)
            .collect();
        let flux = vec![1.0; n_points];
        let flux_err = vec![0.001; n_points];
        LightCurve {
            filename: "flat.csv".to_string(),
            time,
            flux,
            flux_err,
        }
    }

    #[test]
    fn test_bls_recovers_injected_period() {
        // Inject transit with period=3.0 days, depth=0.01
        let lc = make_transit_lc(5000, 3.0, 0.0, 0.03, 0.01, 30.0);
        let periods = generate_periods(0.5, 10.0, 5000);

        let result = bls_search(&lc, &periods, 200, 0.01, 0.05).unwrap();

        // Should recover period within 2%
        let period_error = (result.period - 3.0).abs() / 3.0;
        assert!(
            period_error < 0.02,
            "Period error {:.4} too large (got {:.4}, expected 3.0)",
            period_error,
            result.period
        );
        assert!(result.depth > 0.0, "Depth should be positive");
        assert!(result.power > 0.0, "Power should be positive");
    }

    #[test]
    fn test_bls_recovers_different_period() {
        // Inject transit with period=5.5 days
        let lc = make_transit_lc(8000, 5.5, 0.2, 0.03, 0.015, 55.0);
        let periods = generate_periods(0.5, 15.0, 8000);

        let result = bls_search(&lc, &periods, 200, 0.01, 0.05).unwrap();

        let period_error = (result.period - 5.5).abs() / 5.5;
        assert!(
            period_error < 0.02,
            "Period error {:.4} too large (got {:.4}, expected 5.5)",
            period_error,
            result.period
        );
    }

    #[test]
    fn test_bls_returns_none_for_constant_flux() {
        let lc = make_flat_lc(1000, 30.0);
        let periods = generate_periods(0.5, 10.0, 1000);
        let result = bls_search(&lc, &periods, 200, 0.01, 0.05);
        assert!(result.is_none(), "Constant flux should return None");
    }

    #[test]
    fn test_bls_returns_none_for_too_few_points() {
        let lc = LightCurve {
            filename: "short.csv".to_string(),
            time: vec![1.0, 2.0, 3.0],
            flux: vec![1.0, 0.99, 1.0],
            flux_err: vec![0.001, 0.001, 0.001],
        };
        let periods = generate_periods(0.5, 10.0, 100);
        let result = bls_search(&lc, &periods, 200, 0.01, 0.05);
        assert!(result.is_none(), "Too few points should return None");
    }

    #[test]
    fn test_bls_skips_long_periods() {
        // Time span is 10 days, period > 5 should be skipped
        let lc = make_transit_lc(1000, 2.0, 0.0, 0.03, 0.01, 10.0);
        let periods = vec![8.0, 9.0, 10.0]; // All > time_span/2
        let result = bls_search(&lc, &periods, 200, 0.01, 0.05);
        assert!(result.is_none());
    }

    #[test]
    fn test_estimate_snr_with_transit() {
        let lc = make_transit_lc(5000, 3.0, 0.0, 0.03, 0.01, 30.0);
        let (snr, n_transits) = estimate_snr(&lc, 3.0, 0.0, 0.03);
        assert!(snr > 0.0, "SNR should be positive for transit signal");
        assert!(n_transits >= 2, "Should detect multiple transits over 30 days at P=3d");
    }

    #[test]
    fn test_estimate_snr_empty_light_curve() {
        let lc = LightCurve {
            filename: "empty.csv".to_string(),
            time: vec![],
            flux: vec![],
            flux_err: vec![],
        };
        let (snr, n_transits) = estimate_snr(&lc, 3.0, 0.0, 0.03);
        assert_eq!(snr, 0.0);
        assert_eq!(n_transits, 0);
    }

    #[test]
    fn test_estimate_snr_flat_has_low_snr() {
        let lc = make_flat_lc(1000, 30.0);
        let (snr, _) = estimate_snr(&lc, 3.0, 0.0, 0.03);
        // Flat light curve: depth ≈ 0, so SNR should be near 0
        assert!(snr.abs() < 1.0, "Flat LC should have near-zero SNR, got {}", snr);
    }

    #[test]
    fn test_estimate_snr_phase_wrapping() {
        // Transit at phase 0.97 with duration 0.06 wraps to phase 0.03
        let lc = make_transit_lc(5000, 3.0, 0.97, 0.06, 0.01, 30.0);
        let (snr, _) = estimate_snr(&lc, 3.0, 0.97, 0.06);
        assert!(snr > 0.0, "Phase-wrapped transit should still be detected");
    }

    #[test]
    fn test_generate_periods() {
        let periods = generate_periods(0.5, 20.0, 1000);
        assert_eq!(periods.len(), 1000);
        assert!((periods[0] - 0.5).abs() < 1e-10);
        assert!((periods[999] - 20.0).abs() < 1e-6);
        // Should be monotonically increasing
        for i in 1..periods.len() {
            assert!(periods[i] > periods[i - 1]);
        }
    }

    #[test]
    fn test_compute_phases() {
        let time = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let phases = compute_phases(&time, 2.0, 0.0);
        assert_eq!(phases.len(), 6);
        assert!((phases[0] - 0.0).abs() < 1e-10);
        assert!((phases[1] - 0.5).abs() < 1e-10);
        assert!((phases[2] - 0.0).abs() < 1e-10); // wraps
        assert!((phases[3] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_compute_phases_with_epoch() {
        let time = vec![10.0, 11.0, 12.0];
        let phases = compute_phases(&time, 2.0, 10.0);
        assert!((phases[0] - 0.0).abs() < 1e-10);
        assert!((phases[1] - 0.5).abs() < 1e-10);
        assert!((phases[2] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_median_odd() {
        assert!((median(&[3.0, 1.0, 2.0]) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_median_even() {
        assert!((median(&[4.0, 1.0, 3.0, 2.0]) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_median_empty() {
        assert_eq!(median(&[]), 0.0);
    }

    #[test]
    fn test_median_single() {
        assert!((median(&[42.0]) - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_deeper_transit_has_higher_power() {
        let lc_shallow = make_transit_lc(5000, 3.0, 0.0, 0.03, 0.005, 30.0);
        let lc_deep = make_transit_lc(5000, 3.0, 0.0, 0.03, 0.02, 30.0);
        let periods = generate_periods(0.5, 10.0, 3000);

        let r1 = bls_search(&lc_shallow, &periods, 200, 0.01, 0.05).unwrap();
        let r2 = bls_search(&lc_deep, &periods, 200, 0.01, 0.05).unwrap();

        assert!(
            r2.power > r1.power,
            "Deeper transit should have higher BLS power"
        );
    }
}
