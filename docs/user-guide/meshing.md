# Meshing

MRI2FE provides powerful mesh generation capabilities for creating high-quality finite element meshes from medical imaging data. The meshing module leverages the Computational Geometry Algorithms Library (CGAL) for robust 3D mesh generation.

## Overview

The meshing workflow consists of two main steps:

1. **Format Conversion**: Convert NIfTI images to INR format for CGAL processing
2. **Mesh Generation**: Generate tetrahedral meshes using CGAL's 3D mesh generation algorithms

## Core Functions

### `mesh_from_nifti()`

The primary function for generating meshes from NIfTI files.

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

#### Parameters

- **`filepath`** (str): Path to the input NIfTI file (.nii or .nii.gz)
- **`optimize`** (bool, optional): Enable mesh optimization. Default: `False`
- **`facetAngle`** (float, optional): Minimum angle between facets in degrees. Default: `30.0`
- **`facetSize`** (float, optional): Maximum facet size. Default: `1.0`
- **`facetDistance`** (float, optional): Maximum distance from facet to surface. Default: `4.0`
- **`cellRadiusEdgeRatio`** (float, optional): Maximum radius-edge ratio for tetrahedra. Default: `3.0`
- **`cellSize`** (float, optional): Maximum cell size. Default: `1.0`

#### Returns

- **`meshio.Mesh`**: A meshio mesh object containing nodes and tetrahedral elements

### `nifti_to_inr()`

Converts NIfTI files to INR format for CGAL processing.

```python
from MRI2FE.generate_mesh import nifti_to_inr

inr_path = nifti_to_inr("brain.nii")
```

#### Parameters

- **`filepath`** (str): Path to the input NIfTI file

#### Returns

- **`str`**: Path to the temporary INR file

## Mesh Quality Parameters

### Facet Parameters

- **`facetAngle`**: Controls the smoothness of the surface mesh

  - Lower values (20-30°) produce smoother surfaces but slower generation
  - Higher values (40-50°) produce coarser surfaces but faster generation
  - Recommended range: 20-40°

- **`facetSize`**: Controls the maximum size of surface triangles

  - Smaller values produce finer surface meshes
  - Larger values produce coarser surface meshes
  - Should be scaled with your image resolution

- **`facetDistance`**: Controls how closely the mesh follows the surface
  - Smaller values produce more accurate surface representation
  - Larger values allow more deviation from the surface
  - Recommended: 2-5 times the image voxel size

### Volume Parameters

- **`cellSize`**: Controls the maximum size of tetrahedral elements

  - Smaller values produce finer volume meshes
  - Larger values produce coarser volume meshes
  - Affects computational cost of FE analysis

- **`cellRadiusEdgeRatio`**: Controls the quality of tetrahedral elements
  - Lower values produce higher quality elements
  - Higher values allow lower quality elements
  - Recommended range: 2.0-4.0

## Mesh Optimization

The `optimize` parameter enables CGAL's mesh optimization algorithms:

- **Edge collapse**: Removes unnecessary edges
- **Vertex relocation**: Improves element quality
- **Face removal**: Removes low-quality faces

Optimization improves mesh quality but increases generation time.

## Examples

### Basic Mesh Generation

```python
from MRI2FE.generate_mesh import mesh_from_nifti

# Generate a basic mesh
mesh = mesh_from_nifti("brain.nii")

print(f"Mesh contains {len(mesh.points)} nodes")
print(f"Mesh contains {len(mesh.cells[0][1])} tetrahedral elements")
```

### High-Quality Mesh

```python
# Generate a high-quality mesh for accurate analysis
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=20.0,    # Smooth surfaces
    facetSize=0.5,      # Fine surface mesh
    facetDistance=2.0,  # Accurate surface representation
    cellRadiusEdgeRatio=2.0,  # High-quality elements
    cellSize=0.5        # Fine volume mesh
)
```

### Coarse Mesh for Testing

```python
# Generate a coarse mesh for quick testing
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=False,
    facetAngle=45.0,    # Coarse surfaces
    facetSize=2.0,      # Coarse surface mesh
    facetDistance=8.0,  # Allow surface deviation
    cellRadiusEdgeRatio=4.0,  # Lower quality elements
    cellSize=2.0        # Coarse volume mesh
)
```

### Adaptive Mesh Based on Image Resolution

```python
import ants
from MRI2FE.generate_mesh import mesh_from_nifti

# Load image to get spacing information
img = ants.image_read("brain.nii")
spacing = img.spacing
min_spacing = min(spacing)

# Scale mesh parameters based on image resolution
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=30.0,
    facetSize=min_spacing * 2,      # 2x minimum voxel size
    facetDistance=min_spacing * 4,  # 4x minimum voxel size
    cellRadiusEdgeRatio=3.0,
    cellSize=min_spacing * 2        # 2x minimum voxel size
)
```

## Working with Generated Meshes

### Accessing Mesh Data

```python
# Get node coordinates
nodes = mesh.points
print(f"Node coordinates shape: {nodes.shape}")

# Get element connectivity
elements = mesh.cells[0][1]  # First cell type (tetrahedra)
print(f"Element connectivity shape: {elements.shape}")

# Get mesh bounds
bounds = mesh.bounds
print(f"Mesh bounds: {bounds}")
```

### Converting to FEModel

```python
from MRI2FE.models.femodel import model_from_meshio

# Convert meshio mesh to FEModel
femodel = model_from_meshio(
    mesh=mesh,
    title="Brain Model",
    source="MRI2FE generated"
)

print(f"FE Model has {femodel.metadata['num_nodes']} nodes")
print(f"FE Model has {femodel.metadata['num_elements']} elements")
```

### Saving Meshes

```python
# Save as meshio format
meshio.write("brain_mesh.vtk", mesh)

# Save as LS-DYNA format
femodel.write_lsdyna("brain_model.k")
```

## Troubleshooting

### Common Issues

1. **Memory errors during mesh generation:**

   ```python
   # Reduce mesh density
   mesh = mesh_from_nifti(
       filepath="brain.nii",
       optimize=False,  # Disable optimization
       cellSize=3.0,    # Increase cell size
       facetSize=3.0    # Increase facet size
   )
   ```

2. **Poor mesh quality:**

   ```python
   # Improve mesh quality
   mesh = mesh_from_nifti(
       filepath="brain.nii",
       optimize=True,
       facetAngle=25.0,  # Reduce facet angle
       cellRadiusEdgeRatio=2.5  # Reduce radius-edge ratio
   )
   ```

3. **Long generation times:**
   ```python
   # Speed up generation
   mesh = mesh_from_nifti(
       filepath="brain.nii",
       optimize=False,  # Disable optimization
       facetAngle=40.0,  # Increase facet angle
       cellSize=2.0      # Increase cell size
   )
   ```

### Performance Tips

1. **Start with coarse meshes** for testing and validation
2. **Use appropriate parameters** based on your image resolution
3. **Enable optimization** only for final production meshes
4. **Monitor memory usage** for large datasets

## Advanced Usage

### Batch Processing

```python
import os
from MRI2FE.generate_mesh import mesh_from_nifti

def batch_mesh_generation(input_dir, output_dir):
    """Generate meshes for all NIfTI files in a directory."""

    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith(('.nii', '.nii.gz')):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, f"{filename[:-4]}_mesh.vtk")

            try:
                print(f"Processing {filename}...")
                mesh = mesh_from_nifti(input_path, optimize=True)
                meshio.write(output_path, mesh)
                print(f"Saved mesh to {output_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

# Usage
batch_mesh_generation("input_data/", "output_meshes/")
```

### Custom Mesh Parameters

```python
def create_adaptive_mesh(nifti_path, quality_level="medium"):
    """Create mesh with quality-based parameters."""

    if quality_level == "low":
        params = {
            "optimize": False,
            "facetAngle": 45.0,
            "facetSize": 2.0,
            "facetDistance": 8.0,
            "cellRadiusEdgeRatio": 4.0,
            "cellSize": 2.0
        }
    elif quality_level == "medium":
        params = {
            "optimize": True,
            "facetAngle": 30.0,
            "facetSize": 1.0,
            "facetDistance": 4.0,
            "cellRadiusEdgeRatio": 3.0,
            "cellSize": 1.0
        }
    elif quality_level == "high":
        params = {
            "optimize": True,
            "facetAngle": 20.0,
            "facetSize": 0.5,
            "facetDistance": 2.0,
            "cellRadiusEdgeRatio": 2.0,
            "cellSize": 0.5
        }
    else:
        raise ValueError("quality_level must be 'low', 'medium', or 'high'")

    return mesh_from_nifti(nifti_path, **params)
```

## References

- [CGAL 3D Mesh Generation](https://doc.cgal.org/latest/Mesh_3/)
- [meshio Documentation](https://github.com/nschloe/meshio)
- [NIfTI Format Specification](https://nifti.nimh.nih.gov/nifti-1/)
