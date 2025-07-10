This package contains wrapper functions around Computational Geometry Algorithms Library (CGAL)[@cgalmain] functions for generating a multi-region tetrahedral mesh from a segmented 3D image.

# Example

```python
from MRI2FE import mesh_from_nifti
import meshio


#create mesh
img_path = /path/to/segmented/image.nii

nifti_mesh:meshio.Mesh = mesh_from_nifti(
    img_path,
    optimize=True,
    facetAngle=30,
    facetSize=1.0,
    facetDistance=4.0,
    cellRadiusEdgeRatio=3.0,
    cellSize=1.0
)

#save mesh to abaqus .inp file
nifti_mesh.write("path/to/output",file_format="abaqus")

```

# Reference

::: MRI2FE.mesh_from_nifti