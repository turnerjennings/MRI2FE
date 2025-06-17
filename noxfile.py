import nox
import time
import os

nox.options.reuse_venv = "yes"  # or "yes"

@nox.session(name="test")
def tests(session):
    start_time = time.time()
    session.install("pytest")
    session.install(".")

    #check if test data exist, create if not
    if not all([
        os.path.exists(os.path.join("test", "test_data", "test_stiffness.nii")),
        os.path.exists(os.path.join("test", "test_data", "test_damping_ratio.nii")),
        os.path.exists(os.path.join("test", "test_data", "test_mesh.k"))
    ]):
        print("Generating test data")
        session.run("python","test/create_test_data.py")
    
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
