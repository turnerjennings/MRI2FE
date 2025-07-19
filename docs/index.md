# MRI2FE

Patient-specific Finite Element Model Generation from medical imaging data

## Introduction

This package provides workflows for the generation of finite element (FE) models of the human head from patient-specific magnetic resonance imaging (MRI) and magnetic resonance elastography data. The objective of this toolkit is to provide fast and accurate generation of finite element head models (FEHMs) for a variety of physics and multi-physics analyses. This package leverages industry standard tools including Advanced Normalization Tools (ANTs) for MRI transformation and segmentation [@tustison2021antsx;@KLEINRegistration], and the Computational Geometry Algorithms Library (CGAL) for meshing [@cgalmain;@cgal3d].

## Quick Start

Get up and running quickly with MRI2FE:

```python
import MRI2FE
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import model_from_meshio

# Generate mesh from MRI data
mesh = mesh_from_nifti("brain.nii", optimize=True)

# Create finite element model
femodel = model_from_meshio(mesh, title="Brain Model")
femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])

# Save for simulation
femodel.write_lsdyna("brain_model.k")
```

For a complete workflow example, see the [Quick Start Guide](getting-started/quick-start.md).

## Key Features

### 🧠 **Patient-Specific Modeling**

- Generate FE models directly from MRI data
- Support for various medical imaging formats (NIfTI, DICOM)
- Automatic mesh generation with quality control

### 🔬 **MRE Integration**

- Magnetic Resonance Elastography (MRE) data processing
- Patient-specific material property mapping
- Prony series parameter calculation for viscoelastic materials

### 🏗️ **Advanced Meshing**

- High-quality tetrahedral mesh generation using CGAL
- Configurable mesh quality parameters
- Automatic mesh optimization

### 📊 **Comprehensive I/O**

- Multiple output formats (LS-DYNA, VTK, STL)
- Control keyword generation for simulations
- Post-processing and analysis tools

### 🔧 **Robust Processing**

- Image registration and coregistration
- Error handling and validation
- Batch processing capabilities

## Documentation Structure

### Getting Started

- **[Installation](getting-started/installation.md)** - Install MRI2FE on your system
- **[Quick Start](getting-started/quick-start.md)** - Your first FE model in minutes

### User Guide

- **[User Guide](user-guide.md)** - Complete guide to meshing, model creation, and I/O operations

### Examples

- **[Examples](examples.md)** - Complete workflow examples including basic and MRE integration

## Installation

### Prerequisites

- Python 3.9, 3.10, or 3.11
- CMake 3.15+
- C++ compiler (GCC, Clang, or MSVC)

### Quick Install

**macOS:**

```bash
git clone https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE
./install_mac.sh
```

**Linux:**

```bash
git clone https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE
./install_mac_linux.sh
```

**Windows:**

```cmd
git clone https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE
install_windows.bat
```

For detailed installation instructions, see the [Installation Guide](getting-started/installation.md).

## Basic Usage

### 1. Generate Mesh from MRI

```python
from MRI2FE.generate_mesh import mesh_from_nifti

# Generate high-quality mesh
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=30.0,
    cellSize=1.0
)
```

### 2. Create Finite Element Model

```python
from MRI2FE.models.femodel import model_from_meshio

# Convert mesh to FE model
femodel = model_from_meshio(
    mesh=mesh,
    title="Brain Model",
    source="MRI2FE"
)

# Add material properties
femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])
```

### 3. Export for Simulation

```python
# Save as LS-DYNA format
femodel.write_lsdyna("brain_model.k")
```

## Advanced Features

### MRE Material Property Mapping

```python
from MRI2FE.Pipelines.new_model import NewModel

pipeline = NewModel()

# Register MRE to MRI
registered_data = pipeline.mre_mri_registration(
    mri_geometry_path="brain_mri.nii",
    mre_type="stiffness_damping",
    stiffness_path="stiffness.nii",
    damping_ratio_path="damping.nii"
)

# Calculate material constants
material_constants = pipeline.calculate_material_constants(
    registered_images=registered_data,
    mre_type="stiffness_damping",
    mre_frequency=50.0
)
```

### Multi-Part Models

```python
# Create model with multiple tissue types
femodel = FEModel(title="Multi-Tissue Model")

# Add different parts
femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])
femodel.add_part(2, "skull", [2.0, 0.3, 0.05])
femodel.add_part(3, "cerebrospinal_fluid", [0.1, 0.1, 0.01])
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](development/contributing.md) for details on:

- Setting up a development environment
- Code style and standards
- Testing procedures
- Submitting pull requests

## Citation

If you use MRI2FE in your research, please cite:

```bibtex
@software{mri2fe,
  title={MRI2FE: Patient-specific Finite Element Model Generation},
  author={Jennings, Turner},
  year={2024},
  url={https://github.com/turnerjennings/MRI2FE}
}
```

## Support

- **Documentation**: This site provides comprehensive documentation
- **Issues**: Report bugs and request features on [GitHub](https://github.com/turnerjennings/MRI2FE/issues)
- **Examples**: See the [Examples](examples.md) section for complete workflows

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/turnerjennings/MRI2FE/blob/main/LICENSE) file for details.
