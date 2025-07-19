# API Reference

This page provides API documentation for all MRI2FE modules.

## Table of Contents

- [Core Modules](#core-modules)
- [MRE Processing](#mre-processing)
- [Utilities](#utilities)
- [Output Formats](#output-formats)

## Core Modules

### FEModel

::: MRI2FE.models.femodel.FEModel
handler: python
selection:
members: - **init** - add_nodes - add_elements - add_part - update_centroids - get_part_info - write_lsdyna - **repr**

::: MRI2FE.models.femodel.model_from_meshio
handler: python

### Mesh Generation

::: MRI2FE.generate_mesh.mesh_from_nifti
handler: python

::: MRI2FE.generate_mesh.nifti_to_inr
handler: python

### Utilities

::: MRI2FE.utilities.COM_align
handler: python

::: MRI2FE.utilities.point_cloud_spacing
handler: python

::: MRI2FE.utilities.ants_affine
handler: python

::: MRI2FE.utilities.spatial_map
handler: python

::: MRI2FE.utilities.element_centroids
handler: python

## MRE Processing

### Coregistration

::: MRI2FE.MRE.MRE_coregistration.coregister_MRE_images
handler: python

::: MRI2FE.MRE.MRE_coregistration.segment_MRE_regions
handler: python

::: MRI2FE.MRE.MRE_coregistration.meshio_to_femodel
handler: python

### Material Property Calculation

::: MRI2FE.MRE.calculate_prony.calculate_prony
handler: python

### Mapping

::: MRI2FE.MRE.MRE_mapping.map_MRE_to_mesh
handler: python

### Pipeline

::: MRI2FE.Pipelines.new_model.NewModel
handler: python
selection:
members: - mre_mri_registration - segment_MRE_images - calculate_material_constants - calculate_element_centroids - assign_material_properties

## Utilities

### Image Processing

::: MRI2FE.utilities.COM_align
handler: python

::: MRI2FE.utilities.ants_affine
handler: python

::: MRI2FE.utilities.spatial_map
handler: python

### Geometry and Spatial Analysis

::: MRI2FE.utilities.point_cloud_spacing
handler: python

::: MRI2FE.utilities.element_centroids
handler: python

### Validation Functions

::: MRI2FE.utilities.check_xyz
handler: python

## Output Formats

### LS-DYNA Operations

::: MRI2FE.output.k_file_operations.parse_k_file
handler: python

::: MRI2FE.output.k_file_operations.write_head_k_file
handler: python

::: MRI2FE.output.k_file_operations.edit_control_keyword
handler: python

### Postprocessing

::: MRI2FE.Postprocess.d3_to_nifti.d3_to_displacement
handler: python

::: MRI2FE.Postprocess.d3_to_nifti.grid_to_nifti
handler: python

::: MRI2FE.Postprocess.d3_to_nifti.cloud_to_grid
handler: python

::: MRI2FE.Postprocess.calculate_strain.calculate_strain
handler: python
