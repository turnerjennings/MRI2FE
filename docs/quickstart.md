# Quick Start

Models can be generated quickly from medical imaging data using the built in helper functions.  At it's core, this library contains functions which perform 4 operations:

1. Generating a tetrahedral mesh from a segmented, labeled MRI image.
2. Segmented MRE data into similar regions and calculating prony series viscoelastic parameters for each region.
3. Transforming the MRE data into the segmented MRI space and mapping the segmented MRE material assignments onto the FE mesh.
4. Organizing and writing the data from steps 1-3 to a compatible file format.


## Examples

All steps can be completed at once with a short build script:

```python
import MRI2FE

#define structural MRI paths
labeled_geom_path = "path/to/labeled_image.nii"
geom_roi_mask_path = "path/to/roi_mask.nii"

#define MRE geometry paths
MRE_geometry_paths = ["path/to/30Hzgeom.nii",
                      "path/to/50Hzgeom.nii",
                      "path/to/70Hzgeom.nii"]

MRE_mask_path = "/path/to/MRE_mask.nii"

#define a list of tuples containing the MRE data
#either stiffness/damping ratio or G'/G"
MRE_properties_paths = [
    ("path/to/30Hzstiffness.nii","path/to/30Hzdamping.nii"),
    ("path/to/50Hzstiffness.nii","path/to/50Hzdamping.nii"),
    ("path/to/70Hzstiffness.nii","path/to/70Hzdamping.nii")
]

#model builder workflow: returns model object and writes to output

mdl = (
    MRI2FE.FEModelbuilder()
    .mesh(img_path = labeled_geom_path,
          img_labels = ["region1","region2","region3"])
    .map_mre(target_label = 1,
             MRE_type = "stiffness_damping",
             MRE_geom = MRE_geometry_paths,
             MRE_mask = MRE_mask_path,
             MRE_frequency = [30,50,70],
             MRE_to_transform = MRE_properties_paths)
    .write("/output/path/example.k")
    .build()
)

```

The FEModelbuilder pipeline provides a convenient way to combine all steps of the model generation process.  If additional control over each step is necessary, each step of the build process can be run separately with additional options:

```python
import MRI2FE
import meshio

#create a blank model
mdl = MRI2FE.FEModel(title = "example",
                     source = "example")

#create a mesh with mesh quality constraints
labeled_geom_path = "path/to/geometry.nii"

mesh = MRI2FE.mesh_from_nifti(labeled_geom_path,
                              optimize=True,
                              facetAngle=40,
                              cellRadiusEdgeRatio=2.0)

mdl.from_meshio(mesh)

#load a mesh generated using an external program
external_mesh = meshio.read("path/to/mesh.stl")

mdl.from_meshio(external_mesh)

#transform MRE data using an affine transform

MRE_geometry_paths = ["path/to/30Hzgeom.nii",
                      "path/to/50Hzgeom.nii",
                      "path/to/70Hzgeom.nii"]

MRE_mask_path = "/path/to/MRE_mask.nii"

MRE_properties_paths = [
    ("path/to/30Hzstiffness.nii","path/to/30Hzdamping.nii"),
    ("path/to/50Hzstiffness.nii","path/to/50Hzdamping.nii"),
    ("path/to/70Hzstiffness.nii","path/to/70Hzdamping.nii")
]

_, transformed = MRI2FE.MRE.coregister_MRE_images(
    labeled_geom_path,
    target_label = 4,
    MRE_geom = MRE_geometry_paths,
    type_of_transform = "Affine",
    MRE_to_transform = MRE_properties_paths,
    MRE_frequency = [30,50,70])

#map MRE onto mesh with custom prony series values
mdl = MRI2FE.MRE.map_MRE_to_mesh(
    mdl,
    region_properties = prony_series_values,
    label_background_id = 0,
    region_prefix = "brain"
)

#write output
mdl.write_lsdyna("/save/path/example.k")

```