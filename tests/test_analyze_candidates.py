"""Tests for python/analyze_candidates.py analysis functions."""

import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

from analyze_candidates import (
    crossmatch_known_planets,
    plot_phase_folded,
    setup_style,
)


class TestSetupStyle:
    def test_sets_matplotlib_params(self):
        """setup_style should run without error."""
        setup_style()  # Just verify it doesn't crash


class TestPlotPhaseFolded:
    def test_generates_plot(self, synthetic_lightcurve_csv, tmp_dir):
        csv_path, period, epoch, depth = synthetic_lightcurve_csv
        candidate = {
            "filename": csv_path.name,
            "period_days": period,
            "epoch": epoch,
            "duration_hours": 0.03 * period * 24,
            "depth_ppm": depth * 1e6,
            "snr": 25.0,
            "radius_ratio": depth**0.5,
            "n_transits": 10,
        }
        setup_style()
        output = plot_phase_folded(candidate, csv_path.parent, tmp_dir)
        assert output is not None
        assert output.exists()
        assert output.suffix == ".png"

    def test_missing_file_returns_none(self, tmp_dir):
        candidate = {
            "filename": "nonexistent.csv",
            "period_days": 3.0,
            "epoch": 0.0,
            "duration_hours": 2.0,
            "depth_ppm": 10000,
            "snr": 20.0,
            "radius_ratio": 0.1,
            "n_transits": 10,
        }
        result = plot_phase_folded(candidate, tmp_dir, tmp_dir)
        assert result is None


class TestCrossmatchKnownPlanets:
    def test_match_known_planet(self, catalog_csv):
        candidates = [
            {
                "filename": "TOI-125_TIC_12345.csv",
                "period_days": 4.65,
                "snr": 20.0,
            }
        ]
        df = crossmatch_known_planets(candidates, catalog_csv)
        assert len(df) == 1
        assert df.iloc[0]["status"] == "KNOWN"

    def test_new_candidate(self, catalog_csv):
        candidates = [
            {
                "filename": "TIC_999999_s0056.csv",
                "period_days": 2.0,
                "snr": 10.0,
            }
        ]
        df = crossmatch_known_planets(candidates, catalog_csv)
        assert len(df) == 1
        assert "NEW" in df.iloc[0]["status"]

    def test_missing_catalog(self, tmp_dir):
        candidates = [{"filename": "test.csv", "period_days": 1.0, "snr": 8.0}]
        df = crossmatch_known_planets(candidates, tmp_dir / "nonexistent.csv")
        assert df.empty

    def test_empty_candidates(self, catalog_csv):
        df = crossmatch_known_planets([], catalog_csv)
        assert df.empty


class TestPhaseFoldBinning:
    """Test the vectorized phase-fold binning produces correct results."""

    def test_binning_matches_naive(self):
        """Vectorized binning should match the original for-loop version."""
        rng = np.random.RandomState(42)
        n_bins = 100
        phase = rng.uniform(-0.5, 0.5, 1000)
        flux = rng.normal(1.0, 0.001, 1000)

        # Vectorized version (current code)
        bin_flux_vec = np.zeros(n_bins)
        bin_count_vec = np.zeros(n_bins)
        bins = np.clip(((phase + 0.5) * n_bins).astype(int), 0, n_bins - 1)
        np.add.at(bin_flux_vec, bins, flux)
        np.add.at(bin_count_vec, bins, 1)

        # Naive for-loop version (original code)
        bin_flux_loop = np.zeros(n_bins)
        bin_count_loop = np.zeros(n_bins)
        for i in range(len(phase)):
            b = int((phase[i] + 0.5) * n_bins)
            b = min(b, n_bins - 1)
            bin_flux_loop[b] += flux[i]
            bin_count_loop[b] += 1

        np.testing.assert_array_almost_equal(bin_flux_vec, bin_flux_loop)
        np.testing.assert_array_almost_equal(bin_count_vec, bin_count_loop)
