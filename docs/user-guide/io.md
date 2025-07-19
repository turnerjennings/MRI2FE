# Input/Output

MRI2FE supports various input and output formats for medical imaging data, finite element models, and simulation results. This guide covers all supported formats and how to work with them.

## Overview

MRI2FE handles several data formats:

- **Medical Imaging**: NIfTI (.nii, .nii.gz), INR
- **Finite Element Models**: LS-DYNA (.k), meshio formats
- **Results**: CSV, NIfTI, displacement fields
- **Configuration**: JSON, YAML

## Input Formats

### Medical Imaging Data

#### NIfTI Format

NIfTI (Neuroimaging Informatics Technology Initiative) is the primary format for MRI data.

```python
import ants
from MRI2FE.generate_mesh import mesh_from_nifti

# Load NIfTI file
img = ants.image_read("brain.nii")

# Get image properties
print(f"Image shape: {img.shape}")
print(f"Voxel spacing: {img.spacing}")
print(f"Origin: {img.origin}")

# Generate mesh directly from NIfTI
mesh = mesh_from_nifti("brain.nii")
```

#### INR Format

INR (Inria) format is used internally for CGAL mesh generation.

```python
from MRI2FE.generate_mesh import nifti_to_inr

# Convert NIfTI to INR
inr_path = nifti_to_inr("brain.nii")
print(f"INR file created at: {inr_path}")
```

### Finite Element Models

#### LS-DYNA Format (.k files)

LS-DYNA keyword files are the primary output format for FE models.

```python
from MRI2FE.output.k_file_operations import parse_k_file

# Parse existing LS-DYNA file
element_array, node_array = parse_k_file("existing_model.k")

print(f"Loaded {len(node_array)} nodes")
print(f"Loaded {len(element_array)} elements")
```

#### meshio Formats

meshio supports many formats including VTK, STL, and others.

```python
import meshio
from MRI2FE.models.femodel import model_from_meshio

# Load mesh from various formats
mesh = meshio.read("brain_mesh.vtk")  # VTK format
mesh = meshio.read("brain_mesh.stl")  # STL format
mesh = meshio.read("brain_mesh.msh")  # Gmsh format

# Convert to FEModel
femodel = model_from_meshio(mesh, title="Brain Model")
```

## Output Formats

### Finite Element Models

#### LS-DYNA Format

```python
from MRI2FE.models.femodel import FEModel

# Create model
femodel = FEModel(title="Brain Model")

# Add data...
femodel.add_nodes(node_array=nodes)
femodel.add_elements(element_array=elements)

# Save as LS-DYNA format
femodel.write_lsdyna("output_model.k")
```

#### meshio Formats

```python
import meshio

# Save mesh in various formats
meshio.write("brain_mesh.vtk", mesh)      # VTK format
meshio.write("brain_mesh.stl", mesh)      # STL format
meshio.write("brain_mesh.msh", mesh)      # Gmsh format
meshio.write("brain_mesh.xdmf", mesh)     # XDMF format
```

### Control Keywords

```python
from MRI2FE.output.k_file_operations import edit_control_keyword

# Create control keyword file
edit_control_keyword(
    template="template.k",
    fpath="control.k",
    matprops={
        "brain_tissue": {"Ginf": 1.0, "G1": 0.5, "tau": 0.1},
        "skull": {"Ginf": 2.0, "G1": 0.3, "tau": 0.05}
    },
    includepath="model.k",
    title="Brain Impact Simulation"
)
```

### Results and Analysis

#### Displacement Fields

```python
from MRI2FE.Postprocess.d3_to_nifti import d3_to_displacement

# Convert displacement data to NIfTI
displacement_img = d3_to_displacement(
    displacement_data,
    reference_image,
    output_path="displacement.nii"
)
```

#### Strain Calculations

```python
from MRI2FE.Postprocess.calculate_strain import calculate_strain

# Calculate strain from displacement
strain_tensor = calculate_strain(displacement_field)
```

## Data Processing Functions

### Image Registration

```python
from MRI2FE.MRE.MRE_coregistration import coregister_MRE_images

# Register MRE images to MRI geometry
result = coregister_MRE_images(
    geom="brain_mri.nii",
    gp_list=["storage_modulus.nii"],
    gpp_list=["loss_modulus.nii"],
    imgout="registration_output/"
)
```

### MRE Processing

```python
from MRI2FE.MRE.calculate_prony import calculate_prony

# Calculate Prony series parameters
Ginf, G1, tau = calculate_prony(
    gp=storage_modulus_values,
    gpp=loss_modulus_values,
    w=frequency_array,
    tol=1.0
)
```

### Spatial Mapping

```python
from MRI2FE.utilities import spatial_map

# Create spatial mapping from ANTs image
spatial_data = spatial_map(ants_image)
```

## File Format Specifications

### NIfTI Format

NIfTI files contain:

- **Header**: Image dimensions, spacing, orientation
- **Data**: Voxel intensity values
- **Extensions**: Optional metadata

```python
# NIfTI file properties
img = ants.image_read("brain.nii")
print(f"Dimensions: {img.shape}")
print(f"Spacing: {img.spacing}")
print(f"Origin: {img.origin}")
print(f"Direction: {img.direction}")
print(f"Data type: {img.numpy().dtype}")
```

### LS-DYNA Format

LS-DYNA .k files contain:

- **Keywords**: Control parameters
- **Node data**: Node IDs and coordinates
- **Element data**: Element connectivity and properties
- **Part data**: Material assignments

```python
# LS-DYNA file structure
"""
*KEYWORD
*ELEMENT_SOLID_TET4TOTET10
$#   eid     pid
$#    n1      n2      n3      n4      n5      n6      n7      n8      n9      n10
       1       1       1       2       3       4       0       0       0       0       0
*NODE
#$ nid                x                y               z
       1             0.0             0.0             0.0
       2             1.0             0.0             0.0
*END
"""
```

### INR Format

INR files contain:

- **Header**: Image dimensions, data type, spacing
- **Data**: Raw voxel values

```python
# INR header example
"""
#INRIMAGE-4#{
XDIM=256
YDIM=256
ZDIM=128
VDIM=1
TYPE=float
PIXSIZE=32 bits
CPU=decm
VX=1.0
VY=1.0
VZ=1.0
##}
"""
```

## Examples

### Complete I/O Workflow

```python
import ants
import numpy as np
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import model_from_meshio
from MRI2FE.output.k_file_operations import edit_control_keyword

def complete_workflow(mri_path, output_dir):
    """Complete workflow from MRI to LS-DYNA model."""

    # Step 1: Load MRI data
    print("Loading MRI data...")
    img = ants.image_read(mri_path)
    print(f"Image shape: {img.shape}")

    # Step 2: Generate mesh
    print("Generating mesh...")
    mesh = mesh_from_nifti(mri_path, optimize=True)
    print(f"Mesh: {len(mesh.points)} nodes, {len(mesh.cells[0][1])} elements")

    # Step 3: Create FE model
    print("Creating FE model...")
    femodel = model_from_meshio(mesh, title="Brain Model")
    femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])

    # Step 4: Save model
    print("Saving model...")
    model_path = f"{output_dir}/brain_model.k"
    femodel.write_lsdyna(model_path)

    # Step 5: Create control file
    print("Creating control file...")
    control_path = f"{output_dir}/control.k"
    edit_control_keyword(
        template="template.k",
        fpath=control_path,
        matprops={"brain_tissue": {"Ginf": 1.0, "G1": 0.5, "tau": 0.1}},
        includepath=model_path,
        title="Brain Simulation"
    )

    print("Workflow complete!")
    return femodel

# Usage
model = complete_workflow("brain.nii", "output/")
```

### Batch Processing

```python
import os
import glob
from MRI2FE.generate_mesh import mesh_from_nifti

def batch_process_mri(input_dir, output_dir):
    """Process all NIfTI files in a directory."""

    os.makedirs(output_dir, exist_ok=True)

    # Find all NIfTI files
    nifti_files = glob.glob(os.path.join(input_dir, "*.nii"))
    nifti_files.extend(glob.glob(os.path.join(input_dir, "*.nii.gz")))

    for nifti_file in nifti_files:
        try:
            # Extract filename
            basename = os.path.splitext(os.path.basename(nifti_file))[0]
            if basename.endswith('.nii'):
                basename = os.path.splitext(basename)[0]

            print(f"Processing {basename}...")

            # Generate mesh
            mesh = mesh_from_nifti(nifti_file, optimize=True)

            # Save in multiple formats
            meshio.write(f"{output_dir}/{basename}_mesh.vtk", mesh)
            meshio.write(f"{output_dir}/{basename}_mesh.stl", mesh)

            print(f"Saved mesh for {basename}")

        except Exception as e:
            print(f"Error processing {nifti_file}: {e}")

# Usage
batch_process_mri("input_data/", "output_meshes/")
```

### Data Validation

```python
def validate_nifti_file(filepath):
    """Validate NIfTI file for common issues."""

    try:
        img = ants.image_read(filepath)

        # Check for NaN values
        data = img.numpy()
        if np.any(np.isnan(data)):
            print("Warning: NaN values found in image")

        # Check for infinite values
        if np.any(np.isinf(data)):
            print("Warning: Infinite values found in image")

        # Check image dimensions
        if any(dim <= 0 for dim in img.shape):
            print("Warning: Invalid image dimensions")

        # Check spacing
        if any(sp <= 0 for sp in img.spacing):
            print("Warning: Invalid voxel spacing")

        print(f"Image validation passed: {img.shape}")
        return True

    except Exception as e:
        print(f"Image validation failed: {e}")
        return False

# Usage
validate_nifti_file("brain.nii")
```

## Error Handling

### Common I/O Errors

```python
def robust_file_loading(filepath, file_type="nifti"):
    """Robust file loading with error handling."""

    try:
        if file_type == "nifti":
            return ants.image_read(filepath)
        elif file_type == "mesh":
            return meshio.read(filepath)
        elif file_type == "lsdyna":
            return parse_k_file(filepath)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None
    except PermissionError:
        print(f"Permission denied: {filepath}")
        return None
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

# Usage
img = robust_file_loading("brain.nii", "nifti")
if img is not None:
    print("File loaded successfully")
```

## Performance Considerations

### Large File Handling

```python
def process_large_image(filepath, chunk_size=100):
    """Process large images in chunks."""

    img = ants.image_read(filepath)
    data = img.numpy()

    # Process in chunks to manage memory
    for i in range(0, data.shape[0], chunk_size):
        chunk = data[i:i+chunk_size, :, :]
        # Process chunk...
        print(f"Processed chunk {i//chunk_size + 1}")
```

### Memory-Efficient Mesh Generation

```python
def memory_efficient_mesh(filepath):
    """Generate mesh with memory considerations."""

    # Use conservative parameters for large images
    mesh = mesh_from_nifti(
        filepath=filepath,
        optimize=False,  # Disable optimization to save memory
        cellSize=2.0,    # Larger cells = less memory
        facetSize=2.0    # Larger facets = less memory
    )

    return mesh
```

## Best Practices

1. **Always validate input files** before processing
2. **Use appropriate file formats** for your use case
3. **Handle errors gracefully** with try-catch blocks
4. **Monitor memory usage** for large datasets
5. **Backup important data** before processing
6. **Use descriptive file names** and organize output
7. **Document data provenance** and processing steps

## References

- [NIfTI Format Specification](https://nifti.nimh.nih.gov/nifti-1/)
- [LS-DYNA Keyword Manual](https://www.dynasupport.com/manuals)
- [meshio Documentation](https://github.com/nschloe/meshio)
- [ANTs Documentation](https://antspy.readthedocs.io/)
