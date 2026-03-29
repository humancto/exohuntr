#!/bin/bash
# 🔭 Exoplanet Hunter — Quick Setup
# Run: bash scripts/setup.sh

set -e

echo ""
echo "🔭 Exoplanet Hunter — Setup"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check Rust
if command -v cargo &> /dev/null; then
    echo "✅ Rust $(rustc --version | awk '{print $2}')"
else
    echo "❌ Rust not found. Installing..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
    echo "✅ Rust installed"
fi

# Check Python
if command -v python3 &> /dev/null; then
    echo "✅ Python $(python3 --version | awk '{print $2}')"
else
    echo "❌ Python 3 not found. Please install Python 3.10+"
    exit 1
fi

# Install Python deps
echo ""
echo "📦 Installing Python dependencies..."
pip install lightkurve astroquery pandas numpy matplotlib tqdm --quiet 2>/dev/null || \
pip3 install lightkurve astroquery pandas numpy matplotlib tqdm --quiet

# Build Rust
echo ""
echo "⚙️  Building Rust engine (release mode)..."
cargo build --release

# Create directories
mkdir -p data/lightcurves results/plots

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Setup complete!"
echo ""
echo "🚀 Quick start:"
echo "   make all               # Full pipeline (download + hunt + analyze)"
echo "   make download-candidates  # Download unconfirmed candidates"
echo "   make hunt              # Run BLS detection"
echo "   make analyze           # Generate plots and report"
echo ""
echo "🤖 With Claude Code:"
echo "   claude                 # Launch Claude Code in this directory"
echo "   > 'Hunt for exoplanets in TESS sector 56'"
echo ""
