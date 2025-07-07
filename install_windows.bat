@echo off

setlocal enabledelayedexpansion

set PROJECT_ROOT=%CD%
set VCPKG_TOOLCHAIN=%PROJECT_ROOT%\vcpkg\scripts\buildsystems\vcpkg.cmake

pip install numpy antspyx pybind11

for /f %%i in ('python -c "import pybind11; print(pybind11.get_cmake_dir())"') do set PYBIND11_CMAKE_DIR=%%i

set CMAKE_ARGS=-DCMAKE_TOOLCHAIN_FILE="%VCPKG_TOOLCHAIN%" -DCMAKE_PREFIX_PATH="%PYBIND11_CMAKE_DIR%" -DVCPKG_MANIFEST_MODE=ON -DVCPKG_MANIFEST_DIR="%PROJECT_ROOT%"

python test\create_test_data.py

pip install .