# MRI2FE

Patient-specific Finite Element Model Generation from medical imaging data

## Introduction

This package provides workflows for the generation of finite element (FE) models of the human head from patient-specific magnetic resonance imaging (MRI) and magnetic resonance elastography (MRE) data.  The objective of this toolkit is to provide fast and accurate generation of finite element head models (FEHMs) for a variety of physics and multi-physics analyses.  This package leverages industry standard tools including Advanced Normalization Tools (ANTs) for MRI transformation and segmentation [@tustison2021antsx;@KLEINRegistration], and the Computational Geometry Algorithms Library (CGAL) for meshing [@cgalmain;@cgal3d].

## Installation

Package installation is confirmed working for Windows and MacOS.  To install the package, download or clone the repository to your local machine.  Run the appropriate installation script for your system, which will install all dependencies as well as the package.

```shell
git clone https://github.com/turnerjennings/MRI2FE

cd MRI2FE

#if on windows
./install_windows.bat

#if on mac
./install_mac.sh

#if on linux
./install_linux.sh
```

### Troubleshooting installation

On Windows, after installation you may run into an error stating "DLL load failed".  This issue occurs intermittently on Windows and is due to the dependencies for the CGAL library not linking correctly.  To address the issue, try reinstalling.  If the issue persists, find the location of the DLLs "gmp-10.dll" and "mpfr6.dll" in the ./vcpkg directory and copy them to the installation location.

On Windows, there is an intermittent issue during installation where certain boost libraries cannot be located.  Re-running the installation consistently fixes this issue.