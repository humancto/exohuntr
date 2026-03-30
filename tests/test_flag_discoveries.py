"""Tests for the Phase 2 discovery flagging pipeline."""

import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))
from flag_discoveries import classify_candidate, extract_tic_id


# ===========================================================================
# extract_tic_id
# ===========================================================================

class TestExtractTicId:
    def test_standard_filename(self):
        assert extract_tic_id("TIC_261136679_s0056.csv") == 261136679

    def test_toi_filename(self):
        assert extract_tic_id("TOI_133.01_TIC_219338557.csv") == 219338557

    def test_tic_with_space(self):
        assert extract_tic_id("TIC 261136679_s0056.csv") == 261136679

    def test_no_tic(self):
        assert extract_tic_id("random_star.csv") is None

    def test_lowercase(self):
        assert extract_tic_id("tic_12345_s0001.csv") == 12345

    def test_large_tic_id(self):
        assert extract_tic_id("TIC_999999999_s0070.csv") == 999999999


# ===========================================================================
# classify_candidate
# ===========================================================================

class TestClassifyCandidate:
    def setup_method(self):
        self.toi_tics = {100, 200, 300}
        self.ctoi_tics = {400, 500}
        self.confirmed_tics = {600, 700}

    def test_new_candidate(self):
        status = classify_candidate(
            tic_id=999, planet_score=75, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "NEW"

    def test_known_toi(self):
        status = classify_candidate(
            tic_id=100, planet_score=80, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "KNOWN_TOI"

    def test_known_ctoi(self):
        status = classify_candidate(
            tic_id=400, planet_score=80, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "KNOWN_TOI"

    def test_known_planet(self):
        status = classify_candidate(
            tic_id=600, planet_score=90, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "KNOWN_PLANET"

    def test_confirmed_takes_precedence_over_toi(self):
        """If a TIC is in both confirmed and TOI, it should be KNOWN_PLANET."""
        toi_tics = {600}  # Same as confirmed
        status = classify_candidate(
            tic_id=600, planet_score=90, min_score=60,
            toi_tics=toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "KNOWN_PLANET"

    def test_low_score(self):
        status = classify_candidate(
            tic_id=999, planet_score=40, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "LOW_SCORE"

    def test_no_score_still_new(self):
        """Candidate with no validation score is NEW if not in catalogs."""
        status = classify_candidate(
            tic_id=999, planet_score=None, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "NEW"

    def test_no_tic_id_new_if_high_score(self):
        status = classify_candidate(
            tic_id=None, planet_score=80, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "NEW"

    def test_no_tic_id_low_score(self):
        status = classify_candidate(
            tic_id=None, planet_score=30, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "LOW_SCORE"

    def test_score_at_threshold(self):
        """Score exactly at min_score should pass (>=)."""
        status = classify_candidate(
            tic_id=999, planet_score=60, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "NEW"

    def test_score_just_below_threshold(self):
        status = classify_candidate(
            tic_id=999, planet_score=59, min_score=60,
            toi_tics=self.toi_tics, ctoi_tics=self.ctoi_tics,
            confirmed_tics=self.confirmed_tics,
        )
        assert status == "LOW_SCORE"

    def test_empty_catalogs_all_new(self):
        status = classify_candidate(
            tic_id=100, planet_score=70, min_score=60,
            toi_tics=set(), ctoi_tics=set(), confirmed_tics=set(),
        )
        assert status == "NEW"
