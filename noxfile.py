import nox
import time
import pybind11

nox.options.reuse_venv = "yes"  # or "yes"

@nox.session(name="cpptest")
def cpptest(session):
    import os    
    import shutil
    import sys
    start_time = time.time()

    #define project directories
    project_root = os.getcwd()
    build_dir = "build/test"
    vcpkg_toolchain = os.path.join(
        project_root, "vcpkg/scripts/buildsystems/vcpkg.cmake")

    #delete old build info
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)

    session.install("pybind11")

    # Get the CMake directory from pybind11
    import pybind11
    cmake_dir = pybind11.get_cmake_dir()
    
    # Configure the build
    session.run(
        "cmake",
        "-S", "test",
        "-B", build_dir,
        f"-DCMAKE_TOOLCHAIN_FILE={vcpkg_toolchain}",
        f"-DCMAKE_PREFIX_PATH={cmake_dir}",
        f"-DVCPKG_MANIFEST_MODE=ON",
        f"-DVCPKG_MANIFEST_DIR={project_root}",
        external=True
    )

    # Build the tests
    session.run("cmake", "--build", build_dir, external=True)
    
    # Run the tests
    session.run("ctest", "--test-dir", build_dir, "-C", "Debug", "--output-on-failure", external=True)

    elapsed = time.time() - start_time
    print(f"Session 'cpptest' completed in {elapsed:.2f} seconds")

@nox.session(name="test")
def tests(session):
    start_time = time.time()
    import os    

    project_root = os.path.normpath(os.path.abspath(os.getcwd()))
    vcpkg_toolchain = os.path.normpath(os.path.join(project_root, "vcpkg", "scripts", "buildsystems", "vcpkg.cmake"))

    session.install("pybind11")

    # Get the CMake directory from pybind11
    import pybind11
    cmake_dir = pybind11.get_cmake_dir()

    cmake_args = f'-DCMAKE_TOOLCHAIN_FILE="{vcpkg_toolchain}" '\
        f'-DCMAKE_PREFIX_PATH="{cmake_dir}" ' \
        f'-DVCPKG_MANIFEST_MODE=ON '\
        f'-DVCPKG_MANIFEST_DIR="{project_root}" '
  

    os.environ["CMAKE_ARGS"] = cmake_args

    session.install("scikit-build-core[pyproject]")
    session.install("-v", ".")
    session.run("pytest")

    elapsed = time.time() - start_time
    print(f"Session 'test' completed in {elapsed:.2f} seconds")

@nox.session(name="format")
def format(session):
    start_time = time.time()
    session.install("ruff")
    session.run("ruff","format","src")
    session.run("ruff","format","test")

    elapsed = time.time() - start_time
    print(f"Session 'format' completed in {elapsed:.2f} seconds")

@nox.session(name="lint")
def lint(session):
    start_time = time.time()
    session.install("ruff")
    session.run('ruff','check','src')

    elapsed = time.time() - start_time
    print(f"Session 'lint' completed in {elapsed:.2f} seconds")
