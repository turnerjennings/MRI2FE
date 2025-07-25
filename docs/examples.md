# Examples

This page provides complete examples for common MRI2FE workflows.

## Table of Contents

- [Basic Workflow](#basic-workflow)
- [MRE Integration](#mre-integration)

## Basic Workflow

Complete workflow from MRI data to a finite element model ready for simulation.

### Complete Example

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
        template="template.k",
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

### Parameter Tuning

**High-quality mesh (slower, more accurate):**

```python
mesh_high_quality = mesh_from_nifti(
    filepath="brain.nii",
    optimize=True,
    facetAngle=20.0,    # Smoother surfaces
    facetSize=0.5,      # Finer surface mesh
    facetDistance=2.0,  # More accurate surface
    cellRadiusEdgeRatio=2.0,  # Higher quality elements
    cellSize=0.5        # Finer volume mesh
)
```

**Coarse mesh (faster, less accurate):**

```python
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

### Error Handling

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

        print("Workflow completed successfully!")
        return femodel

    except Exception as e:
        print(f"Workflow failed: {e}")
        return None
```

## MRE Integration

Integrate Magnetic Resonance Elastography (MRE) data to create patient-specific finite element models with accurate material properties.

### Complete MRE Workflow

```python
import numpy as np
import ants
from MRI2FE.Pipelines.new_model import NewModel
from MRI2FE.generate_mesh import mesh_from_nifti
from MRI2FE.models.femodel import model_from_meshio

def create_mre_enhanced_model(
    mri_path,
    stiffness_path,
    damping_ratio_path,
    output_dir,
    mre_frequency=50.0,
    n_segments=5
):
    """
    Create FE model with MRE-derived material properties.

    Args:
        mri_path (str): Path to MRI geometry file
        stiffness_path (str): Path to stiffness map
        damping_ratio_path (str): Path to damping ratio map
        output_dir (str): Output directory
        mre_frequency (float): MRE vibration frequency in Hz
        n_segments (int): Number of MRE segments

    Returns:
        FEModel: Model with MRE-derived material properties
    """

    print("=== MRE-Enhanced Model Creation ===")

    # Initialize pipeline
    pipeline = NewModel()

    # Step 1: Register MRE to MRI
    print("\n1. Registering MRE to MRI geometry...")
    registered_data = pipeline.mre_mri_registration(
        mri_geometry_path=mri_path,
        mre_type="stiffness_damping",
        stiffness_path=stiffness_path,
        damping_ratio_path=damping_ratio_path,
        output_viz_path=f"{output_dir}/registration/"
    )

    # Step 2: Segment MRE regions
    print("\n2. Segmenting MRE regions...")
    segmented_data = pipeline.segment_MRE_images(
        registered_images=registered_data,
        mre_type="stiffness_damping",
        n_segments=n_segments
    )

    # Step 3: Calculate material constants
    print("\n3. Calculating material constants...")
    material_constants = pipeline.calculate_material_constants(
        registered_images=registered_data,
        mre_type="stiffness_damping",
        mre_frequency=mre_frequency
    )

    print("   Material constants calculated:")
    for region, constants in material_constants.items():
        print(f"     {region}: Ginf={constants['Ginf']:.3f}, G1={constants['G1']:.3f}, tau={constants['tau']:.3f}")

    # Step 4: Generate mesh
    print("\n4. Generating mesh...")
    mesh = mesh_from_nifti(
        filepath=mri_path,
        optimize=True,
        facetAngle=30.0,
        facetSize=1.0,
        cellSize=1.0
    )

    # Step 5: Calculate element centroids
    print("\n5. Calculating element centroids...")
    mesh_data = pipeline.calculate_element_centroids("temp_mesh.k")

    # Save mesh temporarily for centroid calculation
    temp_femodel = model_from_meshio(mesh, title="Temp")
    temp_femodel.write_lsdyna("temp_mesh.k")

    # Step 6: Assign material properties
    print("\n6. Assigning material properties...")
    femodel = pipeline.assign_material_properties(
        material_constants=material_constants,
        mesh_data=mesh_data,
        mre_map=registered_data["mu"],  # Use stiffness map
        offset=3,
        output_csv=f"{output_dir}/material_mapping.csv"
    )

    # Step 7: Save final model
    print("\n7. Saving model...")
    model_path = f"{output_dir}/mre_enhanced_model.k"
    femodel.write_lsdyna(model_path)

    # Clean up temporary file
    import os
    if os.path.exists("temp_mesh.k"):
        os.remove("temp_mesh.k")

    print("\n=== MRE Integration Complete ===")
    return femodel

# Usage example
if __name__ == "__main__":
    import os

    # Create output directory
    output_dir = "output/mre_integration"
    os.makedirs(output_dir, exist_ok=True)

    # Run MRE workflow
    model = create_mre_enhanced_model(
        mri_path="data/brain_mri.nii",
        stiffness_path="data/stiffness.nii",
        damping_ratio_path="data/damping_ratio.nii",
        output_dir=output_dir,
        mre_frequency=50.0,
        n_segments=5
    )
```

### MRE Data Validation

```python
def validate_mre_data(stiffness_path, damping_path):
    """Validate MRE data quality."""

    stiffness_img = ants.image_read(stiffness_path)
    damping_img = ants.image_read(damping_path)

    # Check data ranges
    stiffness_data = stiffness_img.numpy()
    damping_data = damping_img.numpy()

    print(f"Stiffness range: [{stiffness_data.min():.3f}, {stiffness_data.max():.3f}]")
    print(f"Damping range: [{damping_data.min():.3f}, {damping_data.max():.3f}]")

    # Check for reasonable values
    if stiffness_data.max() > 100:  # kPa
        print("Warning: Very high stiffness values detected")

    if damping_data.max() > 1.0:
        print("Warning: Damping ratio > 1.0 detected")

    # Check for NaN or infinite values
    if np.any(np.isnan(stiffness_data)) or np.any(np.isnan(damping_data)):
        print("Warning: NaN values found in MRE data")

    if np.any(np.isinf(stiffness_data)) or np.any(np.isinf(damping_data)):
        print("Warning: Infinite values found in MRE data")

# Usage
validate_mre_data("stiffness.nii", "damping_ratio.nii")
```

### Complex Shear Modulus Format

```python
# For complex shear modulus format
registered_data = pipeline.mre_mri_registration(
    mri_geometry_path="brain_mri.nii",
    mre_type="complex_shear",
    complex_modulus_path="complex_modulus.nii"
)

material_constants = pipeline.calculate_material_constants(
    registered_images=registered_data,
    mre_type="complex_shear",
    mre_frequency=50.0
)
```

### Material Property Analysis

```python
def analyze_regional_properties(femodel, material_mapping_csv):
    """Analyze material properties by region."""

    import pandas as pd

    # Load material mapping
    mapping_df = pd.read_csv(material_mapping_csv)

    # Group by region
    regional_stats = mapping_df.groupby('region').agg({
        'Ginf': ['mean', 'std', 'min', 'max'],
        'G1': ['mean', 'std', 'min', 'max'],
        'tau': ['mean', 'std', 'min', 'max']
    })

    print("Regional Material Properties:")
    print(regional_stats)

    return regional_stats

# Usage
stats = analyze_regional_properties(model, "material_mapping.csv")
```

### Robust MRE Workflow

```python
def robust_mre_workflow(mri_path, stiffness_path, damping_path, output_dir):
    """Robust MRE workflow with error handling."""

    try:
        pipeline = NewModel()

        # Validate input files
        for path in [mri_path, stiffness_path, damping_path]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"File not found: {path}")

        # Step 1: Registration with fallback
        try:
            registered_data = pipeline.mre_mri_registration(
                mri_geometry_path=mri_path,
                mre_type="stiffness_damping",
                stiffness_path=stiffness_path,
                damping_ratio_path=damping_path
            )
        except Exception as e:
            print(f"Registration failed: {e}")
            print("Using default material properties...")
            # Fallback to default properties
            return create_basic_model(mri_path, output_dir)

        # Step 2: Material calculation with validation
        try:
            material_constants = pipeline.calculate_material_constants(
                registered_images=registered_data,
                mre_type="stiffness_damping",
                mre_frequency=50.0
            )

            # Validate material constants
            for region, constants in material_constants.items():
                if any(v <= 0 for v in constants.values()):
                    raise ValueError(f"Invalid material constants for region {region}")

        except Exception as e:
            print(f"Material calculation failed: {e}")
            print("Using fallback material properties...")
            material_constants = {"default": {"Ginf": 1.0, "G1": 0.5, "tau": 0.1}}

        # Continue with model creation...
        return create_model_with_materials(mri_path, material_constants, output_dir)

    except Exception as e:
        print(f"MRE workflow failed: {e}")
        return None
```

## Best Practices

1. **Start with basic workflow** before attempting MRE integration
2. **Validate all input data** before processing
3. **Use appropriate mesh parameters** for your use case
4. **Monitor memory usage** for large datasets
5. **Document processing parameters** for reproducibility
6. **Test with small datasets** before processing large files
7. **Backup important data** before processing
8. **Use error handling** for robust workflows
