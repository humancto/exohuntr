"""Shared test fixtures for exohuntr Python tests."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def tmp_dir():
    """Provide a temporary directory that is cleaned up after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def synthetic_lightcurve_csv(tmp_dir):
    """Create a synthetic light curve CSV with a transit signal.

    Returns (csv_path, period, epoch, depth).
    """
    period = 3.0
    epoch = 0.0
    depth = 0.01
    dur_frac = 0.03
    n_points = 5000
    time_span = 30.0

    time = np.linspace(epoch, epoch + time_span, n_points)
    phase = ((time - epoch) / period) % 1.0
    in_transit = (phase < dur_frac / 2) | (phase > 1.0 - dur_frac / 2)
    flux = np.where(in_transit, 1.0 - depth, 1.0)
    # Add small noise
    rng = np.random.RandomState(42)
    flux += rng.normal(0, 0.0001, n_points)
    flux_err = np.full(n_points, 0.001)

    csv_path = tmp_dir / "TOI_999.01_TIC_123456.csv"
    df = pd.DataFrame({"time": time, "flux": flux, "flux_err": flux_err})
    df.to_csv(csv_path, index=False)

    return csv_path, period, epoch, depth


@pytest.fixture
def binary_lightcurve_csv(tmp_dir):
    """Create a synthetic eclipsing binary light curve CSV.

    Has alternating primary/secondary eclipses with different depths
    and a secondary eclipse at phase 0.5.

    Returns (csv_path, period, epoch).
    """
    period = 3.0
    epoch = 0.0
    n_points = 10000
    time_span = 60.0

    time = np.linspace(epoch, epoch + time_span, n_points)
    phase = ((time - epoch) / period) % 1.0
    flux = np.ones(n_points)

    # Primary eclipse at phase 0 (deep)
    primary = (phase < 0.03) | (phase > 0.97)
    flux[primary] = 0.95

    # Secondary eclipse at phase 0.5 (shallower)
    secondary = (phase > 0.47) & (phase < 0.53)
    flux[secondary] = 0.97

    flux_err = np.full(n_points, 0.001)

    csv_path = tmp_dir / "binary_TIC_999999.csv"
    df = pd.DataFrame({"time": time, "flux": flux, "flux_err": flux_err})
    df.to_csv(csv_path, index=False)

    return csv_path, period, epoch


@pytest.fixture
def candidates_json(tmp_dir, synthetic_lightcurve_csv):
    """Create a candidates.json with one candidate matching the synthetic LC."""
    csv_path, period, epoch, depth = synthetic_lightcurve_csv
    candidates = {
        "total_lightcurves": 1,
        "candidates_found": 1,
        "snr_threshold": 6.0,
        "period_range": [0.5, 20.0],
        "candidates": [
            {
                "filename": csv_path.name,
                "period_days": period,
                "epoch": epoch,
                "duration_hours": 0.03 * period * 24,
                "depth_ppm": depth * 1e6,
                "snr": 25.0,
                "n_transits": 10,
                "bls_power": 5.0,
                "radius_ratio": depth ** 0.5,
            }
        ],
    }
    json_path = tmp_dir / "candidates.json"
    with open(json_path, "w") as f:
        json.dump(candidates, f)
    return json_path


@pytest.fixture
def catalog_csv(tmp_dir):
    """Create a mock confirmed exoplanet catalog CSV."""
    df = pd.DataFrame(
        {
            "pl_name": ["TOI-125 b", "WASP-18 b", "HD 21749 b"],
            "hostname": ["TOI-125", "WASP-18", "HD 21749"],
            "pl_orbper": [4.65, 0.94, 35.61],
            "pl_rade": [2.76, 12.4, 2.86],
            "pl_bmasse": [None, 10.3, 22.7],
            "pl_eqt": [None, 2411, 340],
            "st_teff": [5300, 6400, 5200],
            "st_rad": [0.87, 1.23, 0.69],
            "sy_dist": [None, 100.2, 16.3],
            "disc_facility": ["TESS", "SuperWASP", "TESS"],
        }
    )
    path = tmp_dir / "confirmed_exoplanets.csv"
    df.to_csv(path, index=False)
    return path
