# MRE Integration Example

This example demonstrates how to integrate Magnetic Resonance Elastography (MRE) data to create patient-specific finite element models with accurate material properties.

## Overview

MRE integration involves:

1. Registering MRE data to MRI geometry
2. Segmenting MRE regions
3. Calculating material constants (Prony series parameters)
4. Mapping properties to mesh elements
5. Creating a complete FE model

## Complete MRE Workflow

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

    print("   Registration complete")

    # Step 2: Segment MRE regions
    print("\n2. Segmenting MRE regions...")
    segmented_data = pipeline.segment_MRE_images(
        registered_images=registered_data,
        mre_type="stiffness_damping",
        n_segments=n_segments
    )

    print(f"   Created {n_segments} MRE segments")

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

    print(f"   Mesh: {len(mesh.points)} nodes, {len(mesh.cells[0][1])} elements")

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

    print("   Material properties assigned")

    # Step 7: Save final model
    print("\n7. Saving model...")
    model_path = f"{output_dir}/mre_enhanced_model.k"
    femodel.write_lsdyna(model_path)

    print(f"   Model saved to: {model_path}")

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

## Step-by-Step Breakdown

### Step 1: MRE Registration

```python
from MRI2FE.Pipelines.new_model import NewModel

pipeline = NewModel()

# Register MRE data to MRI geometry
registered_data = pipeline.mre_mri_registration(
    mri_geometry_path="brain_mri.nii",
    mre_type="stiffness_damping",
    stiffness_path="stiffness.nii",
    damping_ratio_path="damping_ratio.nii",
    output_viz_path="registration_output/"
)

# Access registered images
stiffness_registered = registered_data["mu"]
damping_registered = registered_data["xi"]
```

### Step 2: MRE Segmentation

```python
# Segment MRE data into regions
segmented_data = pipeline.segment_MRE_images(
    registered_images=registered_data,
    mre_type="stiffness_damping",
    n_segments=5
)

# Access segmentation results
segmentation_map = segmented_data["segmentation"]
prony_parameters = segmented_data["prony_parameters"]
```

### Step 3: Material Constant Calculation

```python
# Calculate Prony series parameters
material_constants = pipeline.calculate_material_constants(
    registered_images=registered_data,
    mre_type="stiffness_damping",
    mre_frequency=50.0
)

# Access material constants
for region, constants in material_constants.items():
    Ginf = constants["Ginf"]  # Long-term shear modulus
    G1 = constants["G1"]      # Short-term shear modulus
    tau = constants["tau"]    # Relaxation time
    print(f"Region {region}: Ginf={Ginf:.3f}, G1={G1:.3f}, tau={tau:.3f}")
```

### Step 4: Property Assignment

```python
# Calculate element centroids for mapping
mesh_data = pipeline.calculate_element_centroids("mesh.k")

# Assign material properties to elements
femodel = pipeline.assign_material_properties(
    material_constants=material_constants,
    mesh_data=mesh_data,
    mre_map=registered_data["mu"],
    offset=3,
    output_csv="material_mapping.csv"
)
```

## Advanced MRE Processing

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

### Custom Segmentation

```python
from MRI2FE.MRE.MRE_coregistration import segment_MRE_regions

# Custom segmentation with specific parameters
segmented_data = segment_MRE_regions(
    SS_img=registered_data["mu"],      # Storage modulus
    DR_img=registered_data["xi"],      # Damping ratio
    n_segs=7                          # Custom number of segments
)
```

### Quality Control

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

## Material Property Analysis

### Regional Analysis

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

### Property Distribution

```python
import matplotlib.pyplot as plt

def plot_property_distribution(material_mapping_csv):
    """Plot distribution of material properties."""

    import pandas as pd

    mapping_df = pd.read_csv(material_mapping_csv)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot Ginf distribution
    axes[0].hist(mapping_df['Ginf'], bins=50, alpha=0.7)
    axes[0].set_title('Long-term Shear Modulus (Ginf)')
    axes[0].set_xlabel('Ginf (kPa)')
    axes[0].set_ylabel('Frequency')

    # Plot G1 distribution
    axes[1].hist(mapping_df['G1'], bins=50, alpha=0.7)
    axes[1].set_title('Short-term Shear Modulus (G1)')
    axes[1].set_xlabel('G1 (kPa)')
    axes[1].set_ylabel('Frequency')

    # Plot tau distribution
    axes[2].hist(mapping_df['tau'], bins=50, alpha=0.7)
    axes[2].set_title('Relaxation Time (tau)')
    axes[2].set_xlabel('tau (s)')
    axes[2].set_ylabel('Frequency')

    plt.tight_layout()
    plt.savefig('material_property_distribution.png', dpi=300)
    plt.show()

# Usage
plot_property_distribution("material_mapping.csv")
```

## Error Handling

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

1. **Validate MRE data** before processing
2. **Use appropriate segmentation** for your tissue types
3. **Check material constants** for physical reasonableness
4. **Document MRE parameters** (frequency, acquisition details)
5. **Compare with literature values** for validation
6. **Use multiple frequency data** when available
7. **Consider tissue anisotropy** for advanced models

## Next Steps

After MRE integration:

1. **Add boundary conditions** for your specific simulation
2. **Validate against experimental data**
3. **Perform sensitivity analysis** on material parameters
4. **Consider multi-scale modeling** approaches
5. **Integrate with other imaging modalities** (DTI, fMRI)

For more advanced workflows, see the [Basic Workflow Example](basic-workflow.md) and [API Reference](../api/mre.md).
