# User Guide

This guide covers all aspects of using MRI2FE for finite element model generation from medical imaging data.

## Table of Contents

- [Meshing](#meshing)
- [Model Creation](#model-creation)
- [Input/Output](#inputoutput)

## Meshing

MRI2FE provides powerful mesh generation capabilities using CGAL for robust 3D mesh generation.

### Core Functions

#### `mesh_from_nifti()`

Generate meshes from NIfTI files:

```python
from MRI2FE.generate_mesh import mesh_from_nifti

mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=30.0,
    facetSize=1.0,
    facetDistance=4.0,
    cellRadiusEdgeRatio=3.0,
    cellSize=1.0
)
```

**Parameters:**

- `filepath`: Path to NIfTI file
- `optimize`: Enable mesh optimization (default: False)
- `facetAngle`: Minimum angle between facets in degrees (default: 30.0)
- `facetSize`: Maximum facet size (default: 1.0)
- `facetDistance`: Maximum distance from facet to surface (default: 4.0)
- `cellRadiusEdgeRatio`: Maximum radius-edge ratio for tetrahedra (default: 3.0)
- `cellSize`: Maximum cell size (default: 1.0)

### Mesh Quality Parameters

- **`facetAngle`**: Controls surface smoothness (20-40° recommended)
- **`facetSize`**: Controls surface mesh density (scale with image resolution)
- **`facetDistance`**: Controls surface accuracy (2-5x voxel size recommended)
- **`cellSize`**: Controls volume mesh density (affects computational cost)
- **`cellRadiusEdgeRatio`**: Controls element quality (2.0-4.0 recommended)

### Examples

**High-quality mesh:**

```python
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=20.0,
    facetSize=0.5,
    cellSize=0.5
)
```

**Coarse mesh for testing:**

```python
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=False,
    facetAngle=45.0,
    facetSize=2.0,
    cellSize=2.0
)
```

## Model Creation

The `FEModel` class handles nodes, elements, parts, materials, and sections.

### Creating Models

#### From meshio Mesh

```python
from MRI2FE.models.femodel import model_from_meshio

femodel = model_from_meshio(
    mesh=mesh,
    title="Brain Model",
    source="MRI2FE"
)
```

#### From Arrays

```python
import numpy as np
from MRI2FE.models.femodel import FEModel

nodes = np.array([
    [1, 0.0, 0.0, 0.0],
    [2, 1.0, 0.0, 0.0],
    [3, 0.0, 1.0, 0.0],
    [4, 0.0, 0.0, 1.0]
])

elements = np.array([
    [1, 1, 1, 2, 3, 4]
])

femodel = FEModel(
    title="Simple Model",
    source="Manual creation",
    nodes=nodes,
    elements=elements
)
```

### Data Structures

**Node Table**: `[node_id, x, y, z]`
**Element Table**: `[element_id, part_id, node1, node2, node3, node4, ...]`

### Adding Parts and Materials

```python
# Add part with material constants
femodel.add_part(
    part_id=1,
    name="brain_tissue",
    material_constants=[1.0, 0.5, 0.1]  # Ginf, G1, tau
)

# Add multiple parts
femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])
femodel.add_part(2, "skull", [2.0, 0.3, 0.05])
```

### Working with MRE Data

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

# Assign properties to mesh
femodel = pipeline.assign_material_properties(
    material_constants=material_constants,
    mesh_data=mesh_data,
    mre_map=mre_image,
    offset=3
)
```

## Input/Output

MRI2FE supports various formats for medical imaging and finite element models.

### Input Formats

#### Medical Imaging

**NIfTI Format:**

```python
import ants
from MRI2FE.generate_mesh import mesh_from_nifti

# Load NIfTI file
img = ants.image_read("brain.nii")
print(f"Image shape: {img.shape}")
print(f"Voxel spacing: {img.spacing}")

# Generate mesh directly
mesh = mesh_from_nifti("brain.nii")
```

#### Finite Element Models

**LS-DYNA Format:**

```python
from MRI2FE.output.k_file_operations import parse_k_file

# Parse existing LS-DYNA file
element_array, node_array = parse_k_file("existing_model.k")
```

**meshio Formats:**

```python
import meshio
from MRI2FE.models.femodel import model_from_meshio

# Load various formats
mesh = meshio.read("brain_mesh.vtk")  # VTK
mesh = meshio.read("brain_mesh.stl")  # STL
mesh = meshio.read("brain_mesh.msh")  # Gmsh

# Convert to FEModel
femodel = model_from_meshio(mesh, title="Brain Model")
```

### Output Formats

#### Finite Element Models

**LS-DYNA Format:**

```python
# Save as LS-DYNA format
femodel.write_lsdyna("output_model.k")
```

**meshio Formats:**

```python
import meshio

# Save in various formats
meshio.write("brain_mesh.vtk", mesh)      # VTK
meshio.write("brain_mesh.stl", mesh)      # STL
meshio.write("brain_mesh.msh", mesh)      # Gmsh
meshio.write("brain_mesh.xdmf", mesh)     # XDMF
```

#### Control Keywords

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

### Data Processing

#### Image Registration

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

#### MRE Processing

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

### Complete Workflow Example

```python
import ants
import numpy as np
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import model_from_meshio
from MRI2FE.output.k_file_operations import edit_control_keyword

def complete_workflow(mri_path, output_dir):
    """Complete workflow from MRI to LS-DYNA model."""

    # Step 1: Load MRI data
    img = ants.image_read(mri_path)
    print(f"Image shape: {img.shape}")

    # Step 2: Generate mesh
    mesh = mesh_from_nifti(mri_path, optimize=True)
    print(f"Mesh: {len(mesh.points)} nodes, {len(mesh.cells[0][1])} elements")

    # Step 3: Create FE model
    femodel = model_from_meshio(mesh, title="Brain Model")
    femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])

    # Step 4: Save model
    model_path = f"{output_dir}/brain_model.k"
    femodel.write_lsdyna(model_path)

    # Step 5: Create control file
    control_path = f"{output_dir}/control.k"
    edit_control_keyword(
        template="template.k",
        fpath=control_path,
        matprops={"brain_tissue": {"Ginf": 1.0, "G1": 0.5, "tau": 0.1}},
        includepath=model_path,
        title="Brain Simulation"
    )

    return femodel

# Usage
model = complete_workflow("brain.nii", "output/")
```

### Error Handling

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
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None
```

### Best Practices

1. **Validate input files** before processing
2. **Use appropriate file formats** for your use case
3. **Handle errors gracefully** with try-catch blocks
4. **Monitor memory usage** for large datasets
5. **Backup important data** before processing
6. **Use descriptive file names** and organize output
7. **Document data provenance** and processing steps
