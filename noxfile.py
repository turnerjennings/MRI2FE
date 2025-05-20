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
    build_dir = "build/test"

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
        f"-DCMAKE_PREFIX_PATH={cmake_dir}",
        external=True
    )

    # Build the tests
    session.run("cmake", "--build", build_dir, external=True)

    python_dir = os.path.dirname(sys.executable)
    test_exe_dir = os.path.join(os.getcwd(), build_dir, "Debug")
    
    # Ensure the directory exists
    os.makedirs(test_exe_dir, exist_ok=True)
    
    # Copy all DLLs from the Python directory
    session.log(f"Copying Python DLLs from {python_dir} to {test_exe_dir}")
    for file in os.listdir(python_dir):
        if file.lower().endswith('.dll'):
            src = os.path.join(python_dir, file)
            dst = os.path.join(test_exe_dir, file)
            try:
                shutil.copy2(src, dst)
                session.log(f"Copied {file}")
            except Exception as e:
                session.log(f"Error copying {file}: {e}")
    
    # Run the tests
    session.run("ctest", "--test-dir", build_dir, "-C", "Debug", "--output-on-failure", external=True)

    elapsed = time.time() - start_time
    print(f"Session 'cpptest' completed in {elapsed:.2f} seconds")

@nox.session(name="test")
def tests(session):
    start_time = time.time()
    session.install("pytest")
    session.install(".")
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
