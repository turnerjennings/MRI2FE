#!/bin/bash
set -e

echo "Installing MRI2FE on macOS..."

# Check if Homebrew is installed
if ! command -v brew &> /dev/null; then
    echo "Homebrew is not installed. Please install Homebrew first:"
    echo "  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
    exit 1
fi

echo "Installing C++ dependencies via Homebrew..."

# Install C++ dependencies
brew install cgal eigen gmp mpfr catch2

echo "Installing Python dependencies..."

# Install Python dependencies
pip install numpy antspyx pybind11 matplotlib lasso-python scipy meshio

echo "Generating test data if needed..."

# Generate test data if needed
python test/create_test_data.py

echo "Installing the MRI2FE package..."

# Install the package
pip install .

echo "Installation complete!" 