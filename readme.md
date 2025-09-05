# MRI2FE

Workflows for automated FE head model generation from MRI/MRE data

## Overview

This package provides tools for the generation of patient-specific finite element (FE) models of the head and brain.  This package currently supports integration of structural data from magnetic resonance imaging (MRI) as well as material data from magnetic resonance elastography (MRE).

## Installation

### Installing from wheel

To install the package, download the appropriate wheel file from the most recent release, and then install on your local machine

```python

pip install ./name-of-wheel.whl

```

### Installing from source

To install the package from source, download or clone the repository to your local machine.  Run the appropriate installation script for your system, which will install all dependencies as well as the package.

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

## Quick Start

All steps can be completed at once with a short build script:

```python
import MRI2FE

#define structural MRI paths
labeled_geom_path = "path/to/labeled_image.nii"
geom_roi_mask_path = "path/to/roi_mask.nii"

#define MRE geometry paths
MRE_geometry_paths = ["path/to/30Hzgeom.nii",
                      "path/to/50Hzgeom.nii",
                      "path/to/70Hzgeom.nii"]

MRE_mask_path = "/path/to/MRE_mask.nii"

#define a list of tuples containing the MRE data
#either stiffness/damping ratio or G'/G"
MRE_properties_paths = [
    ("path/to/30Hzstiffness.nii","path/to/30Hzdamping.nii"),
    ("path/to/50Hzstiffness.nii","path/to/50Hzdamping.nii"),
    ("path/to/70Hzstiffness.nii","path/to/70Hzdamping.nii")
]

#model builder workflow: returns model object and writes to output

mdl = (
    MRI2FE.FEModelbuilder()
    .mesh(img_path = labeled_geom_path,
          img_labels = ["region1","region2","region3"])
    .map_mre(target_label = 1,
             MRE_type = "stiffness_damping",
             MRE_geom = MRE_geometry_paths,
             MRE_mask = MRE_mask_path,
             MRE_frequency = [30,50,70],
             MRE_to_transform = MRE_properties_paths)
    .write("/output/path/example.k")
    .build()
)

```