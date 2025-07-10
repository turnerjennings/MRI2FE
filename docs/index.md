# MRI2FE

Patient-specific Finite Element Model Generation from medical imaging data

## Introduction

This package provides workflows for the generation of finite element (FE) models of the human head from patient-specific magnetic resonance imaging (MRI) and magnetic resonance elastography (MRE) data.  The objective of this toolkit is to provide fast and accurate generation of finite element head models (FEHMs) for a variety of physics and multi-physics analyses.  This package leverages industry standard tools including Advanced Normalization Tools (ANTs) for MRI transformation and segmentation [@tustison2021antsx;@KLEINRegistration], and the Computational Geometry Algorithms Library (CGAL) for meshing [@cgalmain;@cgal3d].

## Installation

Package installation is confirmed working for Windows and MacOS.  To install the package, download or clone the repository to your local machine.  For windows machines, run the script ```install_windows.bat```.  For mac machines, run the ```install_mac.sh``` script.  Each script will install and compile the required dependency libraries before installing the package.

```shell
git clone https://github.com/turnerjennings/MRI2FE

cd MRI2FE

#if on windows
./install_windows.bat

#if on mac
./install_mac.sh
```

## Getting Started

Broadly, this toolkit encompasses functions to perform the following tasks:

1. Generating a tetrahedral FE mesh from a segmented MRI image.
2. Segmenting a MRE image or set of MRE images into discrete regions and calculating prony series material model parameters for each region.
3. Co-registering MRE data to a segmented MRI image and mapping the associated material properties onto a mesh.
4. Writing the resultant mesh to an output file compatible with commercial FE solvers (ABAQUS or LS-DYNA).

The tabs at the top contain instructions on how to use each step in the workflow with your data.  The pipeline tab shows how to use all steps of the process in conjunction to procedurally generate models start-to-finish.  The examples tab contains complete examples of the workflow applied to data. 