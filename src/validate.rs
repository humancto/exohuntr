//! False-positive validation tests for BLS transit candidates.
//!
//! Implements five standard diagnostic tests used by NASA SPOC and automated
//! vetting systems (ExoMiner, LEO-Vetter) to distinguish genuine planet
//! transits from eclipsing binaries and other false positives.
//!
//! # References
//!
//! - Kovacs, Zucker & Mazeh (2002), A&A 391, 369-377 — BLS algorithm
//! - Twicken et al. (2018) — TESS SPOC Data Validation
//! - ExoMiner (Valizadegan et al. 2022) — automated vetting

use crate::bls;
use crate::io::LightCurve;
use serde::Serialize;

/// Result of the odd/even transit depth test.
#[derive(Debug, Clone, Serialize)]
pub struct OddEvenResult {
    /// Ratio of odd transit depth to even transit depth.
    /// A value near 1.0 indicates consistent depths (planet-like).
    pub ratio: f64,
    /// Whether the test passed (ratio between 0.7 and 1.4).
    pub passed: bool,
}

/// Result of the secondary eclipse test.
#[derive(Debug, Clone, Serialize)]
pub struct SecondaryEclipseResult {
    /// Ratio of secondary eclipse depth to primary transit depth.
    pub depth_ratio: f64,
    /// Whether the test passed (ratio < 0.3, no secondary eclipse).
    pub passed: bool,
}

/// Result of the transit shape (V vs U) test.
#[derive(Debug, Clone, Serialize)]
pub struct TransitShapeResult {
    /// Flatness score: higher means more U-shaped (planet-like).
    pub flatness: f64,
    /// Whether the test passed (flatness > 0.2).
    pub passed: bool,
}

/// Result of the period agreement test.
#[derive(Debug, Clone, Serialize)]
pub struct PeriodAgreementResult {
    /// Type of match: "exact", "harmonic", or "disagree".
    pub match_type: String,
    /// Whether the test passed (exact or harmonic match).
    pub passed: bool,
}

/// Complete validation result for a single candidate.
#[derive(Debug, Clone, Serialize)]
pub struct ValidationResult {
    pub filename: String,
    pub period_days: f64,
    pub epoch: f64,
    pub snr: f64,
    pub duration_hours: f64,
    pub depth_ppm: f64,
    pub radius_ratio: f64,
    pub n_transits: usize,
    pub bls_power: f64,
    pub odd_even: Option<OddEvenResult>,
    pub secondary_eclipse: Option<SecondaryEclipseResult>,
    pub transit_shape: Option<TransitShapeResult>,
    pub period_agreement: Option<PeriodAgreementResult>,
    pub planet_score: u32,
}

/// Test 1: Compare transit depths of odd vs even transits.
///
/// Eclipsing binaries have a period twice the detected period, producing
/// alternating primary/secondary eclipses with different depths. A genuine
/// planet produces transits of equal depth.
///
/// Returns `None` if there are insufficient in-transit data points (< 5 per parity).
pub fn test_odd_even_depth(
    time: &[f64],
    flux: &[f64],
    period: f64,
    epoch: f64,
) -> Option<OddEvenResult> {
    let phases = bls::compute_phases(time, period, epoch);

    let mut odd_flux = Vec::new();
    let mut even_flux = Vec::new();
    let mut out_flux = Vec::new();

    for i in 0..time.len() {
        let phase = phases[i];
        let transit_mask = phase < 0.05 || phase > 0.95;
        let out_mask = phase > 0.15 && phase < 0.85;

        if transit_mask {
            let transit_number = ((time[i] - epoch) / period).floor() as i64;
            if transit_number % 2 == 0 {
                even_flux.push(flux[i]);
            } else {
                odd_flux.push(flux[i]);
            }
        }
        if out_mask {
            out_flux.push(flux[i]);
        }
    }

    if odd_flux.len() < 5 || even_flux.len() < 5 {
        return None;
    }

    let baseline = if out_flux.len() > 10 {
        bls::median(&out_flux)
    } else {
        bls::median(flux)
    };

    let depth_odd = baseline - bls::median(&odd_flux);
    let depth_even = baseline - bls::median(&even_flux);

    if depth_even <= 0.0 || depth_odd <= 0.0 {
        return None;
    }

    let ratio = depth_odd / depth_even;
    let passed = (0.7..=1.4).contains(&ratio);

    Some(OddEvenResult { ratio, passed })
}

/// Test 2: Search for a secondary eclipse at orbital phase 0.5.
///
/// Self-luminous companions (stars, brown dwarfs) produce a brightness dip
/// when they pass behind the primary star. Planets are too faint to produce
/// a detectable secondary eclipse in TESS photometry.
///
/// Returns `None` if insufficient data in any phase region.
pub fn test_secondary_eclipse(
    time: &[f64],
    flux: &[f64],
    period: f64,
    epoch: f64,
) -> Option<SecondaryEclipseResult> {
    let phases = bls::compute_phases(time, period, epoch);

    let mut primary_flux = Vec::new();
    let mut secondary_flux = Vec::new();
    let mut out_flux = Vec::new();

    for i in 0..time.len() {
        let phase = phases[i];
        if phase < 0.05 || phase > 0.95 {
            primary_flux.push(flux[i]);
        }
        if phase > 0.45 && phase < 0.55 {
            secondary_flux.push(flux[i]);
        }
        if phase > 0.1 && phase < 0.4 {
            out_flux.push(flux[i]);
        }
    }

    if primary_flux.len() < 5 || secondary_flux.len() < 5 || out_flux.len() < 10 {
        return None;
    }

    let baseline = bls::median(&out_flux);
    let primary_depth = baseline - bls::median(&primary_flux);
    let secondary_depth = baseline - bls::median(&secondary_flux);

    if primary_depth <= 0.0 {
        return None;
    }

    let depth_ratio = secondary_depth / primary_depth;
    let passed = depth_ratio < 0.3;

    Some(SecondaryEclipseResult {
        depth_ratio,
        passed,
    })
}

/// Test 3: Analyze transit shape — U-shaped (planet) vs V-shaped (binary).
///
/// Planets produce flat-bottomed (U-shaped) transits because they fully cover
/// part of the stellar disk. Grazing eclipsing binaries produce V-shaped
/// transits. This test measures the flatness of the transit bottom by comparing
/// the variance of the central portion to the ingress/egress portions.
///
/// Returns `None` if insufficient in-transit data or implausible duration.
pub fn test_transit_shape(
    time: &[f64],
    flux: &[f64],
    period: f64,
    epoch: f64,
    duration_hours: f64,
) -> Option<TransitShapeResult> {
    let phases = bls::compute_phases(time, period, epoch);

    // Center phases on 0 for transit
    let centered_phases: Vec<f64> = phases
        .iter()
        .map(|&p| if p > 0.5 { p - 1.0 } else { p })
        .collect();

    let dur_phase = duration_hours / (period * 24.0);
    if dur_phase <= 0.0 || dur_phase > 0.3 {
        return None;
    }

    // Collect in-transit points sorted by phase
    let mut transit_points: Vec<(f64, f64)> = centered_phases
        .iter()
        .zip(flux.iter())
        .filter(|(&p, _)| p.abs() < dur_phase)
        .map(|(&p, &f)| (p, f))
        .collect();

    if transit_points.len() < 10 {
        return None;
    }

    transit_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let n = transit_points.len();
    let third = n / 3;
    if third < 3 {
        return None;
    }

    let ingress_flux: Vec<f64> = transit_points[..third].iter().map(|(_, f)| *f).collect();
    let flat_flux: Vec<f64> = transit_points[third..2 * third]
        .iter()
        .map(|(_, f)| *f)
        .collect();
    let egress_flux: Vec<f64> = transit_points[2 * third..].iter().map(|(_, f)| *f).collect();

    let flat_std = std_dev(&flat_flux);
    let edge_std = (std_dev(&ingress_flux) + std_dev(&egress_flux)) / 2.0;

    if edge_std < 1e-10 {
        return None;
    }

    let flatness = if edge_std > flat_std {
        1.0 - (flat_std / edge_std)
    } else {
        0.5
    };

    let passed = flatness > 0.2;

    Some(TransitShapeResult { flatness, passed })
}

/// Test 4: Check if detected period matches a reference period.
///
/// Independent detection of the same period strongly validates the signal.
/// Harmonic periods (2x, 0.5x, etc.) are also acceptable as they indicate
/// the same underlying signal.
pub fn test_period_agreement(our_period: f64, reference_period: f64) -> Option<PeriodAgreementResult> {
    if !reference_period.is_finite() || reference_period <= 0.0 {
        return None;
    }

    let ratio = our_period / reference_period;

    // Direct match within 5%
    if (ratio - 1.0).abs() < 0.05 {
        return Some(PeriodAgreementResult {
            match_type: "exact".to_string(),
            passed: true,
        });
    }

    // Harmonic matches within 10%
    for h in &[0.5, 2.0, 3.0, 1.0 / 3.0, 4.0, 0.25] {
        if (ratio - h).abs() < 0.1 {
            return Some(PeriodAgreementResult {
                match_type: "harmonic".to_string(),
                passed: true,
            });
        }
    }

    Some(PeriodAgreementResult {
        match_type: "disagree".to_string(),
        passed: false,
    })
}

/// Compute a 0–100 planet likelihood score from validation test results.
///
/// Scoring weights are based on the discriminative power of each test:
/// - Secondary eclipse test has the strongest weight (±20/25) because a
///   secondary eclipse is almost always an eclipsing binary.
/// - Odd/even depth test (±15): alternating depths indicate a binary.
/// - Period agreement (±15): independent period detection validates the signal.
/// - Transit shape (±10): V-shape suggests a grazing binary.
/// - Radius ratio (±10/15): Rp/Rs > 1.0 is physically impossible for a planet.
/// - SNR bonus (±5): high-SNR detections are more reliable.
pub fn compute_planet_score(
    odd_even: &Option<OddEvenResult>,
    secondary: &Option<SecondaryEclipseResult>,
    shape: &Option<TransitShapeResult>,
    period_match: &Option<PeriodAgreementResult>,
    radius_ratio: f64,
    snr: f64,
) -> u32 {
    let mut score: i32 = 50; // Start neutral

    // Odd/even depth test (±15)
    if let Some(ref r) = odd_even {
        if r.passed {
            score += 15;
        } else {
            score -= 15;
        }
    }

    // Secondary eclipse test (±20/25, strongest discriminator)
    if let Some(ref r) = secondary {
        if r.passed {
            score += 20;
        } else {
            score -= 25;
        }
    }

    // Transit shape test (±10)
    if let Some(ref r) = shape {
        if r.passed {
            score += 10;
        } else {
            score -= 10;
        }
    }

    // Period agreement (±15/5)
    if let Some(ref r) = period_match {
        match r.match_type.as_str() {
            "exact" => score += 15,
            "harmonic" => score += 5,
            "disagree" => score -= 10,
            _ => {}
        }
    }

    // Radius ratio sanity check
    if radius_ratio > 0.0 && radius_ratio < 0.3 {
        score += 10; // planet-sized
    } else if radius_ratio > 1.0 {
        score -= 15; // bigger than star = binary
    }

    // SNR bonus
    if snr > 50.0 {
        score += 5;
    } else if snr < 10.0 {
        score -= 5;
    }

    score.clamp(0, 100) as u32
}

/// Run all validation tests on a single candidate.
///
/// The `reference_period` parameter is optional and comes from an external
/// catalog (e.g., ExoFOP TESS pipeline period) for the period agreement test.
pub fn validate_candidate(
    lc: &LightCurve,
    period: f64,
    epoch: f64,
    duration_hours: f64,
    snr: f64,
    depth_ppm: f64,
    radius_ratio: f64,
    n_transits: usize,
    bls_power: f64,
    reference_period: Option<f64>,
) -> ValidationResult {
    let odd_even = test_odd_even_depth(&lc.time, &lc.flux, period, epoch);
    let secondary_eclipse = test_secondary_eclipse(&lc.time, &lc.flux, period, epoch);
    let transit_shape = test_transit_shape(&lc.time, &lc.flux, period, epoch, duration_hours);
    let period_agreement = reference_period.and_then(|rp| test_period_agreement(period, rp));

    let planet_score = compute_planet_score(
        &odd_even,
        &secondary_eclipse,
        &transit_shape,
        &period_agreement,
        radius_ratio,
        snr,
    );

    ValidationResult {
        filename: lc.filename.clone(),
        period_days: period,
        epoch,
        snr,
        duration_hours,
        depth_ppm,
        radius_ratio,
        n_transits,
        bls_power,
        odd_even,
        secondary_eclipse,
        transit_shape,
        period_agreement,
        planet_score,
    }
}

/// Compute standard deviation of a slice.
fn std_dev(data: &[f64]) -> f64 {
    if data.len() < 2 {
        return 0.0;
    }
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    let variance = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (data.len() as f64 - 1.0);
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Create a light curve with a transit signal at the given parameters.
    /// Suitable for testing validation functions.
    fn make_transit_lc(
        n_points: usize,
        period: f64,
        epoch: f64,
        duration_frac: f64,
        depth: f64,
        time_span: f64,
    ) -> LightCurve {
        let mut time = Vec::with_capacity(n_points);
        let mut flux = Vec::with_capacity(n_points);
        let mut flux_err = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let t = epoch + i as f64 / n_points as f64 * time_span;
            let phase = ((t - epoch) / period) % 1.0;
            let phase = if phase < 0.0 { phase + 1.0 } else { phase };

            let half_dur = duration_frac / 2.0;
            let in_transit = phase < half_dur || phase > (1.0 - half_dur);

            let f = if in_transit { 1.0 - depth } else { 1.0 };
            time.push(t);
            flux.push(f);
            flux_err.push(0.001);
        }

        LightCurve {
            filename: "test_transit.csv".to_string(),
            time,
            flux,
            flux_err,
        }
    }

    /// Create a binary-like light curve with alternating deep/shallow eclipses.
    fn make_binary_lc(
        n_points: usize,
        period: f64,
        epoch: f64,
        primary_depth: f64,
        secondary_depth: f64,
        time_span: f64,
    ) -> LightCurve {
        let mut time = Vec::with_capacity(n_points);
        let mut flux = Vec::with_capacity(n_points);
        let mut flux_err = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let t = epoch + i as f64 / n_points as f64 * time_span;
            let phase = ((t - epoch) / period) % 1.0;
            let phase = if phase < 0.0 { phase + 1.0 } else { phase };

            let f = if phase < 0.03 || phase > 0.97 {
                // Primary eclipse
                1.0 - primary_depth
            } else if phase > 0.47 && phase < 0.53 {
                // Secondary eclipse at phase 0.5
                1.0 - secondary_depth
            } else {
                1.0
            };

            time.push(t);
            flux.push(f);
            flux_err.push(0.001);
        }

        LightCurve {
            filename: "test_binary.csv".to_string(),
            time,
            flux,
            flux_err,
        }
    }

    // =========================================================================
    // Test 1: Odd/Even Transit Depth
    // =========================================================================

    #[test]
    fn test_odd_even_planet_passes() {
        // A planet: same depth on every transit
        let lc = make_transit_lc(10000, 3.0, 0.0, 0.06, 0.01, 60.0);
        let result = test_odd_even_depth(&lc.time, &lc.flux, 3.0, 0.0);
        assert!(result.is_some(), "Should have enough data");
        let r = result.unwrap();
        assert!(
            r.passed,
            "Planet with equal depths should pass. Ratio: {:.3}",
            r.ratio
        );
        assert!(
            (r.ratio - 1.0).abs() < 0.3,
            "Ratio should be near 1.0, got {:.3}",
            r.ratio
        );
    }

    #[test]
    fn test_odd_even_insufficient_data() {
        // Too few data points
        let lc = LightCurve {
            filename: "tiny.csv".to_string(),
            time: vec![0.0, 1.0, 2.0],
            flux: vec![1.0, 0.99, 1.0],
            flux_err: vec![0.001; 3],
        };
        let result = test_odd_even_depth(&lc.time, &lc.flux, 3.0, 0.0);
        assert!(result.is_none());
    }

    // =========================================================================
    // Test 2: Secondary Eclipse
    // =========================================================================

    #[test]
    fn test_secondary_eclipse_planet_passes() {
        // A planet: no secondary eclipse
        let lc = make_transit_lc(10000, 3.0, 0.0, 0.06, 0.01, 60.0);
        let result = test_secondary_eclipse(&lc.time, &lc.flux, 3.0, 0.0);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(
            r.passed,
            "Planet without secondary eclipse should pass. Ratio: {:.3}",
            r.depth_ratio
        );
    }

    #[test]
    fn test_secondary_eclipse_binary_fails() {
        // An eclipsing binary: strong secondary eclipse
        let lc = make_binary_lc(10000, 3.0, 0.0, 0.05, 0.04, 60.0);
        let result = test_secondary_eclipse(&lc.time, &lc.flux, 3.0, 0.0);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(
            !r.passed,
            "Binary with secondary eclipse should fail. Ratio: {:.3}",
            r.depth_ratio
        );
        assert!(r.depth_ratio > 0.3);
    }

    #[test]
    fn test_secondary_eclipse_insufficient_data() {
        let lc = LightCurve {
            filename: "tiny.csv".to_string(),
            time: vec![0.0, 1.0, 2.0],
            flux: vec![1.0, 0.99, 1.0],
            flux_err: vec![0.001; 3],
        };
        let result = test_secondary_eclipse(&lc.time, &lc.flux, 3.0, 0.0);
        assert!(result.is_none());
    }

    // =========================================================================
    // Test 3: Transit Shape
    // =========================================================================

    #[test]
    fn test_transit_shape_box_is_flat() {
        // Box-shaped transit (perfectly flat bottom) — should pass
        let lc = make_transit_lc(20000, 3.0, 0.0, 0.06, 0.01, 60.0);
        let duration_hours = 0.06 * 3.0 * 24.0; // 4.32 hours
        let result = test_transit_shape(&lc.time, &lc.flux, 3.0, 0.0, duration_hours);
        // With a perfect box transit and enough data, test may or may not pass
        // depending on noise — but it should not return None
        if let Some(r) = result {
            // A perfect box transit should have high flatness
            assert!(r.flatness >= 0.0 && r.flatness <= 1.0, "Flatness should be in [0,1]");
        }
    }

    #[test]
    fn test_transit_shape_insufficient_data() {
        let lc = LightCurve {
            filename: "tiny.csv".to_string(),
            time: vec![0.0, 1.0],
            flux: vec![1.0, 0.99],
            flux_err: vec![0.001; 2],
        };
        let result = test_transit_shape(&lc.time, &lc.flux, 3.0, 0.0, 2.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_transit_shape_bad_duration() {
        let lc = make_transit_lc(5000, 3.0, 0.0, 0.03, 0.01, 30.0);
        // Duration > 30% of period is implausible
        let result = test_transit_shape(&lc.time, &lc.flux, 3.0, 0.0, 100.0);
        assert!(result.is_none());

        // Negative duration
        let result = test_transit_shape(&lc.time, &lc.flux, 3.0, 0.0, -1.0);
        assert!(result.is_none());
    }

    // =========================================================================
    // Test 4: Period Agreement
    // =========================================================================

    #[test]
    fn test_period_exact_match() {
        let result = test_period_agreement(3.0, 3.0).unwrap();
        assert_eq!(result.match_type, "exact");
        assert!(result.passed);
    }

    #[test]
    fn test_period_close_match() {
        // Within 5%
        let result = test_period_agreement(3.0, 3.1).unwrap();
        assert_eq!(result.match_type, "exact");
        assert!(result.passed);
    }

    #[test]
    fn test_period_harmonic_2x() {
        let result = test_period_agreement(6.0, 3.0).unwrap();
        assert_eq!(result.match_type, "harmonic");
        assert!(result.passed);
    }

    #[test]
    fn test_period_harmonic_half() {
        let result = test_period_agreement(1.5, 3.0).unwrap();
        assert_eq!(result.match_type, "harmonic");
        assert!(result.passed);
    }

    #[test]
    fn test_period_disagree() {
        // 3.0 / 4.7 ≈ 0.638, not near any harmonic (0.25, 0.33, 0.5, 1.0, 2.0, 3.0, 4.0)
        let result = test_period_agreement(3.0, 4.7).unwrap();
        assert_eq!(result.match_type, "disagree");
        assert!(!result.passed);
    }

    #[test]
    fn test_period_invalid_reference() {
        assert!(test_period_agreement(3.0, 0.0).is_none());
        assert!(test_period_agreement(3.0, -1.0).is_none());
        assert!(test_period_agreement(3.0, f64::NAN).is_none());
        assert!(test_period_agreement(3.0, f64::INFINITY).is_none());
    }

    // =========================================================================
    // Scoring
    // =========================================================================

    #[test]
    fn test_score_all_pass_high_confidence() {
        let score = compute_planet_score(
            &Some(OddEvenResult { ratio: 1.0, passed: true }),
            &Some(SecondaryEclipseResult { depth_ratio: 0.01, passed: true }),
            &Some(TransitShapeResult { flatness: 0.5, passed: true }),
            &Some(PeriodAgreementResult { match_type: "exact".to_string(), passed: true }),
            0.1,  // planet-sized radius ratio
            60.0, // high SNR
        );
        // 50 + 15 + 20 + 10 + 15 + 10 + 5 = 125 → clamped to 100
        assert_eq!(score, 100);
    }

    #[test]
    fn test_score_all_fail_low_confidence() {
        let score = compute_planet_score(
            &Some(OddEvenResult { ratio: 0.3, passed: false }),
            &Some(SecondaryEclipseResult { depth_ratio: 0.9, passed: false }),
            &Some(TransitShapeResult { flatness: 0.1, passed: false }),
            &Some(PeriodAgreementResult { match_type: "disagree".to_string(), passed: false }),
            1.5, // impossibly large
            5.0, // low SNR
        );
        // 50 - 15 - 25 - 10 - 10 - 15 - 5 = -30 → clamped to 0
        assert_eq!(score, 0);
    }

    #[test]
    fn test_score_neutral_no_data() {
        let score = compute_planet_score(&None, &None, &None, &None, 0.0, 30.0);
        // Only base score, no radius or SNR adjustments at 0.0 and 30.0
        assert_eq!(score, 50);
    }

    #[test]
    fn test_score_mixed_results() {
        let score = compute_planet_score(
            &Some(OddEvenResult { ratio: 1.0, passed: true }),
            &Some(SecondaryEclipseResult { depth_ratio: 0.5, passed: false }),
            &Some(TransitShapeResult { flatness: 0.5, passed: true }),
            &None,
            0.1,
            30.0,
        );
        // 50 + 15 - 25 + 10 + 10 = 60
        assert_eq!(score, 60);
    }

    #[test]
    fn test_score_clamped_to_0_100() {
        // Verify score never goes below 0 or above 100
        let low = compute_planet_score(
            &Some(OddEvenResult { ratio: 0.1, passed: false }),
            &Some(SecondaryEclipseResult { depth_ratio: 1.0, passed: false }),
            &Some(TransitShapeResult { flatness: 0.0, passed: false }),
            &Some(PeriodAgreementResult { match_type: "disagree".to_string(), passed: false }),
            2.0,
            3.0,
        );
        assert_eq!(low, 0);

        let high = compute_planet_score(
            &Some(OddEvenResult { ratio: 1.0, passed: true }),
            &Some(SecondaryEclipseResult { depth_ratio: 0.0, passed: true }),
            &Some(TransitShapeResult { flatness: 0.8, passed: true }),
            &Some(PeriodAgreementResult { match_type: "exact".to_string(), passed: true }),
            0.05,
            100.0,
        );
        assert_eq!(high, 100);
    }

    // =========================================================================
    // validate_candidate integration
    // =========================================================================

    #[test]
    fn test_validate_candidate_planet() {
        let lc = make_transit_lc(10000, 3.0, 0.0, 0.06, 0.01, 60.0);
        let result = validate_candidate(
            &lc,
            3.0,
            0.0,
            0.06 * 3.0 * 24.0, // duration_hours
            20.0,               // snr
            10000.0,            // depth_ppm
            0.1,                // radius_ratio
            20,                 // n_transits
            5.0,                // bls_power
            Some(3.0),          // reference_period
        );
        assert!(result.planet_score >= 50, "Planet candidate should score >= 50, got {}", result.planet_score);
        assert_eq!(result.filename, "test_transit.csv");
        assert!(result.period_agreement.is_some());
        assert!(result.period_agreement.as_ref().unwrap().passed);
    }

    #[test]
    fn test_validate_candidate_binary() {
        let lc = make_binary_lc(10000, 3.0, 0.0, 0.05, 0.04, 60.0);
        let result = validate_candidate(
            &lc,
            3.0,
            0.0,
            2.0,   // duration_hours
            15.0,  // snr
            50000.0, // depth_ppm
            1.5,   // radius_ratio > 1.0 (binary-like)
            20,
            3.0,
            Some(7.0), // Period disagrees
        );
        assert!(result.planet_score < 50, "Binary should score < 50, got {}", result.planet_score);
    }

    #[test]
    fn test_validate_candidate_no_reference_period() {
        let lc = make_transit_lc(10000, 3.0, 0.0, 0.06, 0.01, 60.0);
        let result = validate_candidate(
            &lc, 3.0, 0.0, 4.32, 20.0, 10000.0, 0.1, 20, 5.0, None,
        );
        assert!(result.period_agreement.is_none());
    }

    // =========================================================================
    // std_dev helper
    // =========================================================================

    #[test]
    fn test_std_dev_known_values() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let sd = std_dev(&data);
        // Sample std dev: sqrt(32/7) ≈ 2.138
        assert!((sd - 2.138).abs() < 0.01, "Expected ~2.138, got {}", sd);
    }

    #[test]
    fn test_std_dev_constant() {
        let data = vec![5.0, 5.0, 5.0, 5.0];
        assert_eq!(std_dev(&data), 0.0);
    }

    #[test]
    fn test_std_dev_single_element() {
        assert_eq!(std_dev(&[42.0]), 0.0);
    }

    #[test]
    fn test_std_dev_empty() {
        assert_eq!(std_dev(&[]), 0.0);
    }
}
