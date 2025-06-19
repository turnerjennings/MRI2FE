#!/bin/bash
set -e

# Get project root
PROJECT_ROOT=$(pwd)
VCPKG_TOOLCHAIN="$PROJECT_ROOT/vcpkg/scripts/buildsystems/vcpkg.cmake"

# Install Python dependencies
pip install numpy antspyx pybind11

# Get pybind11 CMake directory
PYBIND11_CMAKE_DIR=$(python -c "import pybind11; print(pybind11.get_cmake_dir())")

# Set CMake arguments
export CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=\"$VCPKG_TOOLCHAIN\" -DCMAKE_PREFIX_PATH=\"$PYBIND11_CMAKE_DIR\" -DVCPKG_MANIFEST_MODE=ON -DVCPKG_MANIFEST_DIR=\"$PROJECT_ROOT\""

# Generate test data if needed
python test/create_test_data.py

# Install the package
pip install .