.PHONY: all setup download hunt analyze clean viral test test-rust test-python

# ============================================================================
# 🔭 Exoplanet Hunter — Makefile
# ============================================================================

SECTOR ?= 56
LIMIT ?= 300
SNR ?= 6.0
TOP_N ?= 20
MISSION ?= tess

# Full pipeline: download → hunt → analyze
all: setup download hunt analyze
	@echo ""
	@echo "🏁 Pipeline complete! Check results/REPORT.md"

# Install Python dependencies
setup:
	@echo "📦 Installing dependencies..."
	python3.11 -m pip install -r requirements.txt --quiet

# Download light curves
download:
	@echo "🛰️  Downloading light curves..."
	python3.11 python/download_lightcurves.py \
		--mission $(MISSION) --sector $(SECTOR) \
		--limit $(LIMIT) --catalog

# Download unconfirmed candidates (best for discovery!)
download-candidates:
	@echo "🎯 Downloading unconfirmed TOI candidates..."
	python3.11 python/download_lightcurves.py \
		--candidates-only --limit $(LIMIT) --catalog

# Build and run Rust BLS engine
hunt: target/release/hunt
	@echo "🔭 Running BLS transit detection..."
	./target/release/hunt \
		-i data/lightcurves \
		-o candidates.json \
		--snr-threshold $(SNR) \
		--n-periods 15000

target/release/hunt: src/main.rs Cargo.toml
	@echo "⚙️  Building Rust engine (release mode)..."
	cargo build --release

# Analyze candidates and generate plots
analyze:
	@echo "📊 Analyzing candidates..."
	python3.11 python/analyze_candidates.py \
		--input candidates.json \
		--lightcurves data/lightcurves/ \
		--crossmatch \
		--top-n $(TOP_N)

# Quick hunt: aggressive settings for maximum candidates
aggressive:
	./target/release/hunt \
		-i data/lightcurves \
		-o candidates.json \
		--snr-threshold 4.5 \
		--n-periods 25000 \
		--min-period 0.3 \
		--max-period 40.0

# ============================================================================
# 🧪 Testing
# ============================================================================

# Run all tests (Rust + Python)
test: test-rust test-python

# Run Rust unit tests (66 tests across 4 modules)
test-rust:
	@echo "🧪 Running Rust tests..."
	cargo test

# Run Python tests
test-python:
	@echo "🧪 Running Python tests..."
	python3.11 -m pytest tests/ -v

# Clean generated files
clean:
	rm -rf results/ candidates.json data/lightcurves/
	cargo clean

# ============================================================================
# 🚀 "Make it viral" target
# ============================================================================
viral: all
	@echo ""
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "🚀 VIRAL CHECKLIST:"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	@echo "  1. ✅ Results generated in results/"
	@echo "  2. 📸 Best plots in results/plots/"
	@echo "  3. 📝 Report at results/REPORT.md"
	@echo ""
	@echo "  NEXT:"
	@echo "  • Push to GitHub with a killer README"
	@echo "  • Post best phase-fold plot to r/Astronomy"
	@echo "  • Tweet: 'I built a Rust-powered exoplanet"
	@echo "    hunter and found X candidates in NASA data'"
	@echo "  • Blog post: 'Weekend project: hunting for"
	@echo "    exoplanets with Rust and TESS data'"
	@echo "  • Submit real candidates to ExoFOP/AAVSO"
	@echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
