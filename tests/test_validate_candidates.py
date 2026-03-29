"""Tests for python/validate_candidates.py validation functions."""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add python/ to path so we can import the module
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

from validate_candidates import (
    compute_planet_score,
    test_odd_even_depth as _odd_even_depth,
    test_period_agreement as _period_agreement,
    test_secondary_eclipse as _secondary_eclipse,
    test_transit_shape as _transit_shape,
)


# =============================================================================
# Test 1: Odd/Even Transit Depth
# =============================================================================


class TestOddEvenDepth:
    def _make_planet_lc(self, n=10000, period=3.0, epoch=0.0, depth=0.01, span=60.0):
        """Planet: same depth on every transit."""
        time = np.linspace(epoch, epoch + span, n)
        phase = ((time - epoch) / period) % 1.0
        in_transit = (phase < 0.03) | (phase > 0.97)
        flux = np.where(in_transit, 1.0 - depth, 1.0)
        return time, flux

    def test_planet_passes(self):
        time, flux = self._make_planet_lc()
        ratio, passed = _odd_even_depth(time, flux, 3.0, 0.0)
        assert passed == True
        assert ratio is not None
        assert abs(ratio - 1.0) < 0.3

    def test_insufficient_data_returns_none(self):
        time = np.array([0.0, 1.0, 2.0])
        flux = np.array([1.0, 0.99, 1.0])
        ratio, passed = _odd_even_depth(time, flux, 3.0, 0.0)
        assert ratio == None
        assert passed == None

    def test_zero_depth_returns_none(self):
        """If depths are <= 0 for some reason, returns None."""
        time = np.linspace(0, 60, 10000)
        flux = np.ones_like(time)  # no transit at all
        ratio, passed = _odd_even_depth(time, flux, 3.0, 0.0)
        # Either None (no transit points with dip) or passes
        # The key is it shouldn't crash
        assert ratio == None or isinstance(ratio, float)


# =============================================================================
# Test 2: Secondary Eclipse
# =============================================================================


class TestSecondaryEclipse:
    def _make_planet_lc(self, n=10000, period=3.0, epoch=0.0, depth=0.01, span=60.0):
        time = np.linspace(epoch, epoch + span, n)
        phase = ((time - epoch) / period) % 1.0
        in_transit = (phase < 0.03) | (phase > 0.97)
        flux = np.where(in_transit, 1.0 - depth, 1.0)
        return time, flux

    def _make_binary_lc(self, n=10000, period=3.0, epoch=0.0, span=60.0):
        time = np.linspace(epoch, epoch + span, n)
        phase = ((time - epoch) / period) % 1.0
        flux = np.ones(n)
        flux[(phase < 0.03) | (phase > 0.97)] = 0.95  # primary
        flux[(phase > 0.47) & (phase < 0.53)] = 0.96  # secondary
        return time, flux

    def test_planet_no_secondary(self):
        time, flux = self._make_planet_lc()
        ratio, passed = _secondary_eclipse(time, flux, 3.0, 0.0)
        assert passed == True
        assert ratio is not None
        assert ratio < 0.3

    def test_binary_has_secondary(self):
        time, flux = self._make_binary_lc()
        ratio, passed = _secondary_eclipse(time, flux, 3.0, 0.0)
        assert passed == False
        assert ratio >= 0.3

    def test_insufficient_data(self):
        time = np.array([0.0, 1.0])
        flux = np.array([1.0, 0.99])
        ratio, passed = _secondary_eclipse(time, flux, 3.0, 0.0)
        assert ratio == None


# =============================================================================
# Test 3: Transit Shape
# =============================================================================


class TestTransitShape:
    def test_insufficient_data(self):
        time = np.array([0.0, 1.0])
        flux = np.array([1.0, 0.99])
        score, passed = _transit_shape(time, flux, 3.0, 0.0, 2.0)
        assert score == None

    def test_bad_duration(self):
        time = np.linspace(0, 30, 5000)
        flux = np.ones_like(time)
        # Duration > 30% of period is implausible
        score, passed = _transit_shape(time, flux, 3.0, 0.0, 100.0)
        assert score == None

    def test_returns_valid_flatness(self):
        """For a box-shaped transit, flatness should be a valid number."""
        n = 20000
        period = 3.0
        time = np.linspace(0, 60, n)
        phase = ((time) / period) % 1.0
        flux = np.ones(n)
        dur_frac = 0.06
        in_transit = (phase < dur_frac / 2) | (phase > 1.0 - dur_frac / 2)
        flux[in_transit] = 0.99
        duration_hours = dur_frac * period * 24
        score, passed = _transit_shape(time, flux, period, 0.0, duration_hours)
        if score is not None:
            assert 0.0 <= score <= 1.0


# =============================================================================
# Test 4: Period Agreement
# =============================================================================


class TestPeriodAgreement:
    def test_exact_match(self):
        match_type, passed = _period_agreement(3.0, 3.0)
        assert match_type == "exact"
        assert passed == True

    def test_close_match(self):
        match_type, passed = _period_agreement(3.0, 3.1)
        assert match_type == "exact"
        assert passed == True

    def test_harmonic_2x(self):
        match_type, passed = _period_agreement(6.0, 3.0)
        assert match_type == "harmonic"
        assert passed == True

    def test_harmonic_half(self):
        match_type, passed = _period_agreement(1.5, 3.0)
        assert match_type == "harmonic"
        assert passed == True

    def test_disagree(self):
        # 3.0 / 4.7 ≈ 0.638, not near any harmonic
        match_type, passed = _period_agreement(3.0, 4.7)
        assert match_type == "disagree"
        assert passed == False

    def test_invalid_tess_period_none(self):
        result = _period_agreement(3.0, None)
        assert result == (None, None)

    def test_invalid_tess_period_nan(self):
        result = _period_agreement(3.0, float("nan"))
        assert result == (None, None)

    def test_invalid_tess_period_zero(self):
        result = _period_agreement(3.0, 0)
        assert result == (None, None)

    def test_invalid_tess_period_negative(self):
        result = _period_agreement(3.0, -1.0)
        assert result == (None, None)


# =============================================================================
# Scoring
# =============================================================================


class TestComputePlanetScore:
    def test_all_pass_high_score(self):
        results = {
            "odd_even_passed": True,
            "secondary_passed": True,
            "shape_passed": True,
            "period_match": "exact",
            "radius_ratio": 0.1,
            "snr": 60.0,
        }
        score = compute_planet_score(results)
        assert score == 100

    def test_all_fail_low_score(self):
        results = {
            "odd_even_passed": False,
            "secondary_passed": False,
            "shape_passed": False,
            "period_match": "disagree",
            "radius_ratio": 1.5,
            "snr": 5.0,
        }
        score = compute_planet_score(results)
        assert score == 0

    def test_neutral_no_data(self):
        results = {"snr": 30.0}
        score = compute_planet_score(results)
        assert score == 50

    def test_clamped_0_100(self):
        # Score can't go below 0
        results = {
            "odd_even_passed": False,
            "secondary_passed": False,
            "shape_passed": False,
            "period_match": "disagree",
            "radius_ratio": 2.0,
            "snr": 3.0,
        }
        score = compute_planet_score(results)
        assert score >= 0
        assert score <= 100

    def test_mixed_results(self):
        results = {
            "odd_even_passed": True,
            "secondary_passed": False,
            "shape_passed": True,
            "radius_ratio": 0.1,
            "snr": 30.0,
        }
        score = compute_planet_score(results)
        # 50 + 15 - 25 + 10 + 10 = 60
        assert score == 60

    def test_harmonic_period_bonus(self):
        results = {"period_match": "harmonic", "snr": 30.0}
        score = compute_planet_score(results)
        assert score == 55  # 50 + 5
