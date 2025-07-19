# Basic Workflow Example

This example demonstrates the complete workflow from MRI data to a finite element model ready for simulation.

## Overview

The basic workflow consists of:

1. Loading MRI data
2. Generating a mesh
3. Creating a finite element model
4. Adding material properties
5. Exporting for simulation

## Complete Example

```python
import numpy as np
import ants
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import model_from_meshio
from MRI2FE.output.k_file_operations import edit_control_keyword

def create_brain_model(mri_path, output_dir):
    """
    Create a complete brain FE model from MRI data.

    Args:
        mri_path (str): Path to MRI NIfTI file
        output_dir (str): Directory for output files

    Returns:
        FEModel: The created finite element model
    """

    print("=== MRI2FE Basic Workflow ===")

    # Step 1: Load and validate MRI data
    print("\n1. Loading MRI data...")
    img = ants.image_read(mri_path)
    print(f"   Image shape: {img.shape}")
    print(f"   Voxel spacing: {img.spacing}")
    print(f"   Data range: [{img.numpy().min():.2f}, {img.numpy().max():.2f}]")

    # Step 2: Generate mesh
    print("\n2. Generating mesh...")
    mesh = mesh_from_nifti(
        filepath=mri_path,
        optimize=True,
        facetAngle=30.0,
        facetSize=1.0,
        facetDistance=4.0,
        cellRadiusEdgeRatio=3.0,
        cellSize=1.0
    )

    print(f"   Mesh generated: {len(mesh.points)} nodes, {len(mesh.cells[0][1])} elements")

    # Step 3: Create finite element model
    print("\n3. Creating finite element model...")
    femodel = model_from_meshio(
        mesh=mesh,
        title="Brain FE Model",
        source="MRI2FE Basic Workflow"
    )

    # Step 4: Add material properties
    print("\n4. Adding material properties...")
    femodel.add_part(
        part_id=1,
        name="brain_tissue",
        material_constants=[1.0, 0.5, 0.1]  # Example viscoelastic parameters
    )

    # Step 5: Save model
    print("\n5. Saving model...")
    model_path = f"{output_dir}/brain_model.k"
    femodel.write_lsdyna(model_path)
    print(f"   Model saved to: {model_path}")

    # Step 6: Create control file
    print("\n6. Creating control file...")
    control_path = f"{output_dir}/control.k"
    edit_control_keyword(
        template="template.k",  # You'll need to provide this template
        fpath=control_path,
        matprops={
            "brain_tissue": {
                "Ginf": 1.0,    # Long-term shear modulus
                "G1": 0.5,      # Short-term shear modulus
                "tau": 0.1      # Relaxation time
            }
        },
        includepath=model_path,
        title="Brain Impact Simulation"
    )
    print(f"   Control file saved to: {control_path}")

    print("\n=== Workflow Complete ===")
    print(f"Model statistics:")
    print(f"  - Nodes: {femodel.metadata['num_nodes']}")
    print(f"  - Elements: {femodel.metadata['num_elements']}")
    print(f"  - Parts: {len(femodel.part_info)}")

    return femodel

# Usage example
if __name__ == "__main__":
    import os

    # Create output directory
    output_dir = "output/basic_workflow"
    os.makedirs(output_dir, exist_ok=True)

    # Run workflow
    model = create_brain_model(
        mri_path="data/brain_mri.nii",
        output_dir=output_dir
    )
```

## Step-by-Step Breakdown

### Step 1: Loading MRI Data

```python
import ants

# Load the MRI image
img = ants.image_read("brain.nii")

# Validate the image
print(f"Image dimensions: {img.shape}")
print(f"Voxel spacing: {img.spacing}")
print(f"Data type: {img.numpy().dtype}")

# Check for common issues
data = img.numpy()
if np.any(np.isnan(data)):
    print("Warning: NaN values found in image")
if np.any(np.isinf(data)):
    print("Warning: Infinite values found in image")
```

### Step 2: Mesh Generation

```python
from MRI2FE.generate_mesh import mesh_from_nifti

# Generate mesh with quality parameters
mesh = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,              # Enable mesh optimization
    facetAngle=30.0,           # Smooth surface mesh
    facetSize=1.0,             # Maximum facet size
    facetDistance=4.0,         # Surface accuracy
    cellRadiusEdgeRatio=3.0,   # Element quality
    cellSize=1.0               # Maximum element size
)

# Access mesh data
nodes = mesh.points
elements = mesh.cells[0][1]  # Tetrahedral elements
print(f"Generated {len(nodes)} nodes and {len(elements)} elements")
```

### Step 3: Model Creation

```python
from MRI2FE.models.femodel import model_from_meshio

# Convert mesh to FE model
femodel = model_from_meshio(
    mesh=mesh,
    title="Brain Model",
    source="MRI2FE"
)

# Add part information
femodel.add_part(
    part_id=1,
    name="brain_tissue",
    material_constants=[1.0, 0.5, 0.1]
)
```

### Step 4: Model Export

```python
# Save as LS-DYNA format
femodel.write_lsdyna("brain_model.k")

# Create control file for simulation
from MRI2FE.output.k_file_operations import edit_control_keyword

edit_control_keyword(
    template="template.k",
    fpath="control.k",
    matprops={
        "brain_tissue": {
            "Ginf": 1.0,
            "G1": 0.5,
            "tau": 0.1
        }
    },
    includepath="brain_model.k",
    title="Brain Simulation"
)
```

## Parameter Tuning

### Mesh Quality vs. Performance

```python
# High-quality mesh (slower, more accurate)
mesh_high_quality = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=20.0,    # Smoother surfaces
    facetSize=0.5,      # Finer surface mesh
    facetDistance=2.0,  # More accurate surface
    cellRadiusEdgeRatio=2.0,  # Higher quality elements
    cellSize=0.5        # Finer volume mesh
)

# Coarse mesh (faster, less accurate)
mesh_coarse = mesh_from_nifti(
    filepath="brain.nii",
    optimize=False,
    facetAngle=45.0,    # Coarser surfaces
    facetSize=2.0,      # Coarser surface mesh
    facetDistance=8.0,  # Less accurate surface
    cellRadiusEdgeRatio=4.0,  # Lower quality elements
    cellSize=2.0        # Coarser volume mesh
)
```

### Adaptive Parameters

```python
def adaptive_mesh_parameters(img):
    """Generate mesh parameters based on image properties."""

    spacing = img.spacing
    min_spacing = min(spacing)

    return {
        "optimize": True,
        "facetAngle": 30.0,
        "facetSize": min_spacing * 2,      # 2x minimum voxel size
        "facetDistance": min_spacing * 4,  # 4x minimum voxel size
        "cellRadiusEdgeRatio": 3.0,
        "cellSize": min_spacing * 2        # 2x minimum voxel size
    }

# Use adaptive parameters
img = ants.image_read("brain.nii")
params = adaptive_mesh_parameters(img)
mesh = mesh_from_nifti("brain.nii", **params)
```

## Error Handling

```python
def robust_workflow(mri_path, output_dir):
    """Robust workflow with error handling."""

    try:
        # Step 1: Validate input
        if not os.path.exists(mri_path):
            raise FileNotFoundError(f"MRI file not found: {mri_path}")

        img = ants.image_read(mri_path)
        if img.numpy().size == 0:
            raise ValueError("MRI file is empty")

        # Step 2: Generate mesh with fallback
        try:
            mesh = mesh_from_nifti(mri_path, optimize=True)
        except MemoryError:
            print("Memory error, trying with conservative parameters...")
            mesh = mesh_from_nifti(
                mri_path,
                optimize=False,
                cellSize=3.0,
                facetSize=3.0
            )

        # Step 3: Create model
        femodel = model_from_meshio(mesh, title="Brain Model")
        femodel.add_part(1, "brain_tissue", [1.0, 0.5, 0.1])

        # Step 4: Save with validation
        model_path = f"{output_dir}/brain_model.k"
        femodel.write_lsdyna(model_path)

        if not os.path.exists(model_path):
            raise RuntimeError("Failed to save model file")

        print("Workflow completed successfully!")
        return femodel

    except Exception as e:
        print(f"Workflow failed: {e}")
        return None
```

## Validation and Quality Checks

```python
def validate_model(femodel):
    """Validate the generated model."""

    issues = []

    # Check node data
    if femodel.node_table is None or len(femodel.node_table) == 0:
        issues.append("No nodes in model")

    # Check element data
    if femodel.element_table is None or len(femodel.element_table) == 0:
        issues.append("No elements in model")

    # Check for duplicate node IDs
    if femodel.node_table is not None:
        node_ids = femodel.node_table[:, 0]
        if len(node_ids) != len(set(node_ids)):
            issues.append("Duplicate node IDs found")

    # Check element connectivity
    if femodel.element_table is not None and femodel.node_table is not None:
        node_ids = set(femodel.node_table[:, 0])
        for i, element in enumerate(femodel.element_table):
            node_refs = element[2:]  # Skip element_id and part_id
            for node_ref in node_refs:
                if node_ref not in node_ids:
                    issues.append(f"Element {i+1} references non-existent node {node_ref}")

    return issues

# Usage
model = create_brain_model("brain.nii", "output/")
issues = validate_model(model)
if issues:
    print("Model validation issues:", issues)
else:
    print("Model validation passed!")
```

## Next Steps

After completing the basic workflow:

1. **Add MRE data** for patient-specific material properties
2. **Create multi-part models** for different tissue types
3. **Optimize mesh parameters** for your specific use case
4. **Add boundary conditions** and loading conditions
5. **Run simulations** using LS-DYNA or other FE solvers

See the [MRE Integration Example](mre-integration.md) for advanced workflows with material property mapping.
