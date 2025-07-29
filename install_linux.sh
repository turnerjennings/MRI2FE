#!/bin/bash

set -e  # Exit on any error

#install prereqs
sudo apt update
sudo apt install -y libcgal-dev libeigen3-dev libgmp-dev libmpfr-dev catch2 autoconf automake libtool

PROJECT_ROOT=$(pwd)
VCPKG_TOOLCHAIN="$PROJECT_ROOT/vcpkg/scripts/buildsystems/vcpkg.cmake"

pip install numpy antspyx pybind11

PYBIND11_CMAKE_DIR=$(python -c "import pybind11; print(pybind11.get_cmake_dir())")

export CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$VCPKG_TOOLCHAIN -DCMAKE_PREFIX_PATH=$PYBIND11_CMAKE_DIR -DVCPKG_MANIFEST_MODE=ON -DVCPKG_MANIFEST_DIR=$PROJECT_ROOT"

pip install .

python test/create_test_data.py