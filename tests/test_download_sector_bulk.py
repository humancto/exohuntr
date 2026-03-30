"""Tests for the Phase 2 bulk sector download pipeline."""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))
from download_sector_bulk import fetch_toi_tic_ids


# ===========================================================================
# fetch_toi_tic_ids
# ===========================================================================

class TestFetchToiTicIds:
    @patch("download_sector_bulk.pd.read_csv")
    def test_returns_tic_ids(self, mock_read_csv):
        """Test that TOI and cTOI TIC IDs are merged."""
        # First call = TOI list, second = cTOI list
        toi_df = pd.DataFrame({"TIC ID": [100, 200, 300]})
        ctoi_df = pd.DataFrame({"TIC ID": [300, 400, 500]})
        mock_read_csv.side_effect = [toi_df, ctoi_df]

        result = fetch_toi_tic_ids()

        assert isinstance(result, set)
        assert result == {100, 200, 300, 400, 500}

    @patch("download_sector_bulk.pd.read_csv")
    def test_handles_toi_failure(self, mock_read_csv):
        """If TOI fetch fails, should still try cTOI."""
        ctoi_df = pd.DataFrame({"TIC ID": [400, 500]})
        mock_read_csv.side_effect = [Exception("Network error"), ctoi_df]

        result = fetch_toi_tic_ids()
        assert 400 in result
        assert 500 in result

    @patch("download_sector_bulk.pd.read_csv")
    def test_handles_both_failure(self, mock_read_csv):
        """If both fetches fail, returns empty set."""
        mock_read_csv.side_effect = [Exception("fail"), Exception("fail")]

        result = fetch_toi_tic_ids()
        assert result == set()

    @patch("download_sector_bulk.pd.read_csv")
    def test_deduplicates(self, mock_read_csv):
        """Same TIC ID in both catalogs should appear once."""
        toi_df = pd.DataFrame({"TIC ID": [100, 100, 200]})
        ctoi_df = pd.DataFrame({"TIC ID": [200, 300]})
        mock_read_csv.side_effect = [toi_df, ctoi_df]

        result = fetch_toi_tic_ids()
        assert result == {100, 200, 300}
