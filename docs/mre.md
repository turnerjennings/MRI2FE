# MRE Co-registration

This package contains functionality to convert magnetic resonance elastography (MRE) data into finite element (FE) friendly material formulations, and map the MRE onto an associated FE model.  The functionality can broadly be broken down into four steps:

1. Co-registration of MRE to an associated magnetic resonance imaging (MRI) geometry.  
2. Discretization of MRE into a clustered set of regions to reduce the required number of material constants.
3. Conversion of MRE data represented by complex shear moduli or shear stiffness and damping ratio at different frequencies to a prony series representation compatible with FE prony series or generalized Maxwell material models.
4. Mapping of MRE data to an associated mesh.

# Example

To be inserted.

# API Reference

::: MRI2FE.MRE.coregister_MRE_images

::: MRI2FE.MRE.segment_MRE_regions

::: MRI2FE.MRE.calculate_prony

::: MRI2FE.MRE.map_MRE_to_mesh