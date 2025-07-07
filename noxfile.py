import nox
import time
import os
import pybind11
import platform

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
    is_macos = platform.system() == "Darwin"
    
    if is_macos:
        vcpkg_toolchain = None
        cmake_args = []
    else:
        vcpkg_toolchain = os.path.join(
            project_root, "vcpkg/scripts/buildsystems/vcpkg.cmake")
        cmake_args = [
            f"-DCMAKE_TOOLCHAIN_FILE={vcpkg_toolchain}",
            f"-DVCPKG_MANIFEST_MODE=ON",
            f"-DVCPKG_MANIFEST_DIR={project_root}"
        ]

    # delete old build info
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)

    session.install("pybind11")

    # Get the CMake directory from pybind11
    import pybind11
    cmake_dir = pybind11.get_cmake_dir()
    cmake_args.append(f"-DCMAKE_PREFIX_PATH={cmake_dir}")

    #check if test data already exists, generate if not
    if not os.path.exists(os.path.join(project_root,
                                       "test",
                                       "test_data",
                                       "test_concentric_spheres.nii")):
        print("Test data does not exist, generating...")
        session.install("numpy")
        session.install("antspyx")
        session.run("python", f"{os.path.join(project_root, 'test', 'generate_test_niftis.py')}")


    if not os.path.exists(os.path.join(project_root,
                                       "test",
                                       "test_data",
                                       "test_concentric_spheres.inr")):
        print("Test data does not exist, generating...")
        session.install("numpy")
        session.install("antspyx")
        session.run("python", f"{os.path.join(project_root, 'test', 'generate_test_niftis.py')}")
    
    # Configure the build
    session.run(
        "cmake",
        "-S", "test",
        "-B", build_dir,
        *cmake_args,
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
    import platform

    project_root = os.path.normpath(os.path.abspath(os.getcwd()))
    is_macos = platform.system() == "Darwin"
    
    session.install("pytest")
    session.install("pybind11")

    # check if test data exist, create if not
    if not all([
        os.path.exists(os.path.join("test", "test_data", "test_stiffness.nii")),
        os.path.exists(os.path.join("test", "test_data", "test_damping_ratio.nii")),
        os.path.exists(os.path.join("test", "test_data", "test_mesh.k")),
        os.path.exists(os.path.join("test","test_data","test_concentric_spheres.nii")),
        os.path.exists(os.path.join("test","test_data","test_concentric_spheres.inr"))
    ]):
        print("Generating test data")
        
        session.install("numpy")
        session.install("antspyx")
        session.run("python", f"{os.path.join(project_root, 'test', 'create_test_data.py')}")

    if is_macos:
        session.run("bash", "install_mac.sh")
    else:
        vcpkg_toolchain = os.path.normpath(os.path.join(project_root, "vcpkg", "scripts", "buildsystems", "vcpkg.cmake"))
        
        import pybind11
        cmake_dir = pybind11.get_cmake_dir()

        cmake_args = f'-DCMAKE_TOOLCHAIN_FILE="{vcpkg_toolchain}" '\
            f'-DCMAKE_PREFIX_PATH="{cmake_dir}" ' \
            f'-DVCPKG_MANIFEST_MODE=ON '\
            f'-DVCPKG_MANIFEST_DIR="{project_root}" '
      
        os.environ["CMAKE_ARGS"] = cmake_args
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
