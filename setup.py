#!/usr/bin/env python3
"""Setup script for MRI2FE package."""

import os
import sys
from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Get the long description from the README file
def read_readme():
    with open("README.md", encoding="utf-8") as f:
        return f.read()

# Get version from git or set default
def get_version():
    try:
        import setuptools_scm
        return setuptools_scm.get_version()
    except ImportError:
        return "0.1.0"

# Define C++ extensions
ext_modules = [
    Pybind11Extension(
        "MRI2FE._core",
        ["src/MRI2FE/_core.cpp"],
        include_dirs=[
            "include",
            "src/MRI2FE/include",
        ],
        language="c++",
        cxx_std=17,
    ),
]

# Setup configuration
setup(
    name="MRI2FE",
    version=get_version(),
    description="Finite Element Model Generation from MRI/MRE Data",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/MRI2FE",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/MRI2FE/issues",
        "Source": "https://github.com/yourusername/MRI2FE",
        "Documentation": "https://github.com/yourusername/MRI2FE#readme",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords=["finite-element", "mri", "mre", "biomechanics", "medical-imaging"],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.21.0",
        "antspyx>=0.3.0",
        "pybind11>=2.10.0",
        "matplotlib>=3.5.0",
        "lasso-python>=0.3.0",
        "scipy>=1.7.0",
        "meshio>=5.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "nox>=2022.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
        ],
        "test": [
            "pytest>=6.0",
            "pytest-cov>=3.0.0",
        ],
    },
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    include_package_data=True,
    zip_safe=False,
) 