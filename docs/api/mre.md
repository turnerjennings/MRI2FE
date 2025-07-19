# MRE Processing

This page provides API documentation for Magnetic Resonance Elastography (MRE) processing modules.

## Coregistration

::: MRI2FE.MRE.MRE_coregistration.coregister_MRE_images
handler: python

::: MRI2FE.MRE.MRE_coregistration.segment_MRE_regions
handler: python

::: MRI2FE.MRE.MRE_coregistration.meshio_to_femodel
handler: python

## Material Property Calculation

::: MRI2FE.MRE.calculate_prony.calculate_prony
handler: python

## Mapping

::: MRI2FE.MRE.MRE_mapping.map_MRE_to_mesh
handler: python

## Pipeline

::: MRI2FE.Pipelines.new_model.NewModel
handler: python
selection:
members: - mre_mri_registration - segment_MRE_images - calculate_material_constants - calculate_element_centroids - assign_material_properties
