# Quick Start Guide

This guide will walk you through creating your first finite element model using MRI2FE. We'll cover the basic workflow from MRI data to a complete FE model.

## Prerequisites

Before starting, ensure you have:

- MRI2FE installed (see [Installation](installation.md))
- MRI data in NIfTI format (.nii or .nii.gz)
- Optional: MRE data for material property assignment

## Basic Workflow

### 1. Import Required Modules

```python
import MRI2FE
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import FEModel
from MRI2FE.utilities import COM_align
```

### 2. Generate Mesh from MRI Data

```python
# Load your MRI data and generate a mesh
mesh = mesh_from_nifti(
    filepath="path/to/your/mri_data.nii",
    optimize=True,
    facetAngle=30.0,
    facetSize=1.0,
    cellSize=1.0
)

print(f"Generated mesh with {len(mesh.points)} nodes and {len(mesh.cells[0][1])} elements")
```

### 3. Create Finite Element Model

```python
# Convert mesh to FEModel
femodel = FEModel(
    title="My Brain Model",
    source="MRI2FE generated",
    nodes=mesh.points,
    elements=mesh.cells[0][1]
)

# Add part information
femodel.add_part(part_id=1, name="brain", material_constants=[1.0, 0.5, 0.1])

print(f"Created FE model with {femodel.metadata['num_nodes']} nodes")
```

### 4. Save the Model

```python
# Save as LS-DYNA format
femodel.write_lsdyna("output_model.k")
```

## Complete Example

Here's a complete example that demonstrates the full workflow:

```python
import numpy as np
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import FEModel
from MRI2FE.utilities import COM_align

def create_basic_brain_model(mri_path, output_path):
    """
    Create a basic brain FE model from MRI data.

    Args:
        mri_path (str): Path to MRI NIfTI file
        output_path (str): Path for output LS-DYNA file
    """

    # Step 1: Generate mesh
    print("Generating mesh from MRI data...")
    mesh = mesh_from_nifti(
        filepath=mri_path,
        optimize=True,
        facetAngle=30.0,
        facetSize=1.0,
        facetDistance=4.0,
        cellRadiusEdgeRatio=3.0,
        cellSize=1.0
    )

    # Step 2: Prepare node and element data
    nodes = np.column_stack([
        np.arange(1, len(mesh.points) + 1),  # Node IDs
        mesh.points  # Coordinates
    ])

    elements = np.column_stack([
        np.arange(1, len(mesh.cells[0][1]) + 1),  # Element IDs
        np.ones(len(mesh.cells[0][1])),  # Part ID (all elements in part 1)
        mesh.cells[0][1] + 1  # Node connectivity (1-indexed)
    ])

    # Step 3: Create FE model
    print("Creating finite element model...")
    femodel = FEModel(
        title="Brain FE Model",
        source="MRI2FE",
        nodes=nodes,
        elements=elements
    )

    # Step 4: Add material properties
    femodel.add_part(
        part_id=1,
        name="brain_tissue",
        material_constants=[1.0, 0.5, 0.1]  # Example values
    )

    # Step 5: Save model
    print(f"Saving model to {output_path}...")
    femodel.write_lsdyna(output_path)

    print("Model creation complete!")
    return femodel

# Usage
if __name__ == "__main__":
    model = create_basic_brain_model(
        mri_path="data/brain_mri.nii",
        output_path="output/brain_model.k"
    )
```

## Advanced Features

### Working with MRE Data

If you have MRE data for material property assignment:

```python
from MRI2FE.Pipelines.new_model import NewModel

# Initialize the pipeline
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

### Mesh Optimization

You can control mesh quality with various parameters:

```python
# High-quality mesh (slower generation)
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=20.0,  # Smaller angle = smoother surfaces
    facetSize=0.5,    # Smaller size = finer mesh
    cellSize=0.5      # Smaller size = more elements
)

# Coarse mesh (faster generation)
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=False,
    facetAngle=45.0,
    facetSize=2.0,
    cellSize=2.0
)
```

## Next Steps

Now that you have the basics, explore:

- [Meshing Guide](../user-guide/meshing.md) - Detailed mesh generation options
- [Model Creation](../user-guide/model-creation.md) - Advanced FE model features
- [Input/Output](../user-guide/io.md) - Working with different file formats
- [Examples](../examples/) - Complete workflow examples

## Troubleshooting

### Common Issues

1. **Memory errors during mesh generation:**

   - Reduce `cellSize` and `facetSize` parameters
   - Use `optimize=False` for faster processing

2. **Poor mesh quality:**

   - Adjust `facetAngle` (lower values = smoother surfaces)
   - Enable `optimize=True`
   - Check input image quality

3. **File format issues:**
   - Ensure MRI data is in NIfTI format
   - Check file paths are correct
   - Verify file permissions

For more detailed troubleshooting, see the [Installation Guide](installation.md).
